"""Numerical parity with Bioconductor limma 3.68.5 (reference values from tests/data/make_reference.R)."""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import pytest

import pylimma

DATA = pathlib.Path(__file__).parent / "data"


def _read(name):
    return pd.read_csv(DATA / name, sep="\t", float_precision="round_trip")


def _close(a, b, rtol=1e-8, atol=0.0):
    np.testing.assert_allclose(np.asarray(a, dtype=float), np.asarray(b, dtype=float), rtol=rtol, atol=atol)


@pytest.fixture(scope="module")
def data():
    X = pd.read_csv(DATA / "counts.tsv", sep="\t", index_col=0)
    X = X.iloc[_read("keep.tsv")["keep"].to_numpy()]
    lib = X.columns.str.replace(r"_(KO|NT)$", "", regex=True)
    cond = pd.Categorical(X.columns.str.split("_").str[-1], ["NT", "KO"])
    D = pd.get_dummies(pd.DataFrame({"lib": pd.Categorical(lib), "cond": cond}), drop_first=True).astype(float)
    D.insert(0, "Intercept", 1.0)
    eff = _read("eff_lib.tsv")["lib"].to_numpy()
    return X, cond, D.to_numpy(), eff


@pytest.fixture(scope="module")
def ebayes_fit(data):
    X, _, D, eff = data
    v = pylimma.voom(X, D, lib_size=eff)
    return v, pylimma.eBayes(pylimma.lmFit(v, D))          # lmFit(v, design), as in R


def test_voom(ebayes_fit):
    v, _ = ebayes_fit
    ref = _read("voom.tsv")
    _close(v.E.ravel(order="F"), ref["E"], rtol=1e-12)
    _close(v.weights.ravel(order="F"), ref["W"], rtol=1e-6)


def test_ebayes(ebayes_fit):
    _, f = ebayes_fit
    ref, sc = _read("ebayes.tsv"), _read("ebayes_scalar.tsv")
    k = f.coefficients.shape[1] - 1
    _close(f.df_prior, sc["df_prior"][0], rtol=1e-6)
    _close(f.s2_prior, sc["s2_prior"][0], rtol=1e-6)
    _close(f.coefficients[:, k], ref["coef"], rtol=1e-6, atol=1e-9)
    _close(f.t[:, k], ref["t"], rtol=1e-6, atol=1e-9)
    _close(f.p_value[:, k], ref["p"], rtol=1e-6)
    _close(f.s2_post, ref["s2post"], rtol=1e-6)
    _close(f.var_prior[k], sc["var_prior"][0], rtol=1e-6)
    _close(f.lods[:, k], ref["B"], rtol=1e-6, atol=1e-8)


def test_toptable_sorted_by_b(ebayes_fit):
    _, f = ebayes_fit
    k = f.coefficients.shape[1] - 1
    tt = pylimma.topTable(f, coef=k, number=np.inf)
    ref = _read("toptable.tsv")
    assert list(tt["gene"]) == list(ref["gene"])           # gene names carried from the counts
    _close(tt["B"], ref["B"], rtol=1e-6, atol=1e-8)
    _close(tt["adj.P.Val"], ref["adj"], rtol=1e-6)


def test_continuous_covariate_t_is_not_clipped(data):
    """stdev.coef.lim only bounds the B-statistic's prior; python-limma 0.1.0 clipped it into t."""
    X, cond, _, eff = data
    age = _read("age.tsv")["age"].to_numpy()
    D2 = np.column_stack([np.ones(len(cond)), np.asarray(cond) == "KO", age]).astype(float)
    f = pylimma.eBayes(pylimma.lmFit(pylimma.voom(X, D2, lib_size=eff), D2))
    ref = _read("ebayes_age.tsv")
    assert (ref["su"] < 0.1).all()
    _close(f.t[:, 2], ref["t"], rtol=1e-6, atol=1e-9)
    _close(f.p_value[:, 2], ref["p"], rtol=1e-6)
    _close(f.lods[:, 2], ref["B"], rtol=1e-6, atol=1e-8)


def test_treat(data, ebayes_fit):
    _, _, D, _ = data
    v, _ = ebayes_fit
    tr = pylimma.treat(pylimma.lmFit(v, D), lfc=1)
    ref = _read("treat.tsv")
    k = D.shape[1] - 1
    _close(tr.t[:, k], ref["t"], rtol=1e-6, atol=1e-9)
    _close(tr.p_value[:, k], ref["p"], rtol=1e-6)
    assert tr.lods is None

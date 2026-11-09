"""pylimma run on every benchmark task, mirroring 02_run_limma.R step by step.

Usage: python 03_run_pylimma.py <tasks_dir> <out_dir> <n_timing_reps> [task_glob]
The count preprocessing (filterByExpr, TMM) uses pyedger, which matches edgeR
to ~1e-15 on these data; voom, lmFit, eBayes and treat are pylimma's.
"""
from __future__ import annotations

import copy
import importlib.metadata
import json
import os
import pathlib
import sys
import time

import numpy as np
import pandas as pd

import pyedger
import pylimma


def run_pipeline(X: pd.DataFrame, group: pd.Categorical) -> dict:
    tm = {}
    t0 = time.perf_counter()
    y = pyedger.DGEList(X, group=group)
    keep = np.asarray(pyedger.filterByExpr(y))
    y = pyedger.normLibSizes(y.subset(keep, keep_lib_sizes=False))
    design = np.column_stack([np.ones(len(group)), np.asarray(group) == "Tumor"]).astype(float)
    tm["filter_norm"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    v = pylimma.voom(y.counts, design, lib_size=pyedger.effective_lib_sizes(y))
    tm["voom"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    fit = pylimma.lmFit(v, design)
    efit = pylimma.eBayes(copy.copy(fit))        # eBayes/treat set attributes on the fit they get
    tm["lmFit_eBayes"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    tfit = pylimma.treat(copy.copy(fit), lfc=1)  # reuses the fit, as treat(fit) in R
    tm["treat"] = time.perf_counter() - t0
    return dict(y=y, v=v, keep=keep, efit=efit, tfit=tfit, time=tm)


def bh(p):
    from scipy.stats import false_discovery_control
    return false_discovery_control(np.asarray(p, dtype=float), method="bh")


def engine() -> str:
    return f"threads{os.environ.get('NUMBA_NUM_THREADS', 'all')}"


def main():
    tasks_dir, out_dir, reps = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), int(sys.argv[3])
    pattern = sys.argv[4] if len(sys.argv) > 4 else "*"
    out_dir.mkdir(parents=True, exist_ok=True)
    for td in sorted(tasks_dir.glob(pattern)):
        if not (td / "counts.tsv.gz").exists():
            continue
        od = out_dir / td.name
        if (od / "timing.json").exists():
            continue
        od.mkdir(exist_ok=True)
        X = pd.read_csv(td / "counts.tsv.gz", sep="\t", index_col=0)
        s = pd.read_csv(td / "samples.tsv", sep="\t", index_col=0)
        group = pd.Categorical(s.loc[X.columns, "group"], categories=["Normal", "Tumor"])
        times = []
        for _ in range(reps):
            res = run_pipeline(X, group)
            times.append(res["time"])
        v, e, t = res["v"], res["efit"], res["tfit"]
        lods = getattr(e, "lods", None)
        out = pd.DataFrame({
            "gene": res["y"].counts.index,
            "AveExpr": e.Amean,
            "w_mean": v.weights.mean(axis=1),
            "w_first": v.weights[:, 0],
            "E_first": v.E[:, 0],
            "sigma": e.sigma,
            "s2post": e.s2_post,
            "eb_logFC": e.coefficients[:, 1],
            "eb_t": e.t[:, 1],
            "eb_PValue": e.p_value[:, 1],
            "eb_FDR": bh(e.p_value[:, 1]),
            "eb_B": np.nan if lods is None else np.asarray(lods)[:, 1],
            "tr_t": t.t[:, 1],
            "tr_PValue": t.p_value[:, 1],
            "tr_FDR": bh(t.p_value[:, 1]),
        })
        out.to_csv(od / "results.tsv.gz", sep="\t", index=False, float_format="%.17g")
        f = lambda x: x.tolist() if isinstance(x, np.ndarray) else float(x)
        scal = dict(n_keep=int(res["keep"].sum()), norm_factors=res["y"].samples["norm.factors"].tolist(),
                    df_prior=f(e.df_prior), s2_prior=f(e.s2_prior))
        (od / "scalars.json").write_text(json.dumps(scal))
        (od / "timing.json").write_text(json.dumps(dict(engine=engine(), pylimma=getattr(pylimma, "__version__", "0.1.0"),
                                                        pyedger=pyedger.__version__, reps=times)))
        print(td.name, engine(), f"{sum(times[-1].values()):.2f}s", flush=True)


if __name__ == "__main__":
    main()

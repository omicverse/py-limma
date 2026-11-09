"""Figures, tables and the numbers quoted in the pylimma technical report.

Usage: python 05_figures.py <runs_dir> <fgsea.tsv> <tasks_dir> <out_dir>
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

RUNS, FGSEA, TASKS, OUT = map(pathlib.Path, sys.argv[1:5])
OUT.mkdir(parents=True, exist_ok=True)
ASSETS = pathlib.Path(__file__).resolve().parent.parent / "assets"

C_R, C_PY = "#eb6834", "#2a78d6"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = LinearSegmentedColormap.from_list("seq", ["#fcfcfb", "#cde2fb", "#6da7ec", "#256abf", "#0d366b"])
QCOL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#008300"]
mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42, "savefig.dpi": 300,
})
FDR, LFC = 0.05, 1.0
TESTS = {"eb": "eBayes moderated t", "tr": "treat (lfc = 1)"}


def load(engine, task):
    return pd.read_csv(RUNS / engine / task / "results.tsv.gz", sep="\t")


def calls(df, test):
    sig = df[f"{test}_FDR"] <= FDR
    if test == "eb":
        sig &= df["eb_logFC"].abs() >= LFC
    return np.where(sig & (df["eb_logFC"] > 0), 1, np.where(sig, -1, 0))


def rel_err(a, b, floor=1e-250):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.abs(a - b) / np.maximum(np.abs(a), floor)


def sci(x):
    if x == 0:
        return "0"
    m, e = f"{x:.1e}".split("e")
    return f"${m} \\times 10^{{{int(e)}}}$"


full = sorted(p.name for p in (RUNS / "limma").glob("full_*"))
cohorts = [t.split("_", 1)[1] for t in full]
subs = sorted(p.name for p in (RUNS / "limma").glob("sub*"))
ERRK = {"AveExpr": "AveExpr", "w_mean": "voom weights (gene mean)", "sigma": "residual SD",
        "eb_logFC": "logFC", "eb_PValue": "P (eBayes)", "eb_B": "B-statistic", "tr_PValue": "P (treat)"}

rows, errs = [], {k: [] for k in ERRK}
for task in full + subs:
    r, p = load("limma", task), load("pylimma", task)
    assert (r.gene.to_numpy() == p.gene.to_numpy()).all(), task
    meta = json.load(open(TASKS / task / "meta.json"))
    sr, sp = (json.load(open(RUNS / e / task / "scalars.json")) for e in ("limma", "pylimma"))
    row = dict(task=task, kind=meta["kind"], n=meta.get("n"), samples=meta["n_samples"], genes_kept=len(r),
               dfprior_rel=rel_err([sr["df_prior"]], [sp["df_prior"]])[0],
               s2prior_rel=rel_err([sr["s2_prior"]], [sp["s2_prior"]])[0],
               w_maxrel=rel_err(r.w_first, p.w_first).max())
    for t in TESTS:
        cr, cp = calls(r, t), calls(p, t)
        row[f"{t}_deg_R"], row[f"{t}_discordant"] = int((cr != 0).sum()), int((cr != cp).sum())
        row[f"{t}_p_maxrel"] = np.nanmax(rel_err(r[f"{t}_PValue"], p[f"{t}_PValue"]))
        row[f"{t}_fdr_flip"] = int(((r[f"{t}_FDR"] <= FDR) != (p[f"{t}_FDR"] <= FDR)).sum())
    rows.append(row)
    if meta["kind"] == "full":
        for k in errs:
            errs[k].append(rel_err(r[k], p[k]))
tab = pd.DataFrame(rows)
tab.to_csv(OUT / "comparison_per_task.tsv", sep="\t", index=False)
tim = pd.DataFrame({e: [sum(json.load(open(RUNS / e / t / "timing.json"))["reps"][0].values()) for t in full]
                    for e in ("limma", "pylimma")}, index=cohorts)
tim.insert(0, "samples", tab.set_index("task").loc[full, "samples"].to_numpy())
tim.to_csv(OUT / "timing_full.tsv", sep="\t")
stages = pd.DataFrame({e: json.load(open(RUNS / e / "full_BRCA" / "timing.json"))["reps"][0] for e in ("limma", "pylimma")})
stages.to_csv(OUT / "timing_stages_BRCA.tsv", sep="\t")

# ======================= Figure 1 =======================
fig = plt.figure(figsize=(7.2, 6.9))
gs = fig.add_gridspec(3, 1, height_ratios=[1.3, 1.0, 1.0], hspace=0.38)


def panel(ax, s, dx=-0.02):
    ax.text(dx, 1.03, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right", color=INK)


row1 = gs[0].subgridspec(1, 2, width_ratios=[0.16, 1.0], wspace=0.04)
logo = row1[0].subgridspec(2, 1, height_ratios=[1.25, 0.55], hspace=0.0)
for k, f in enumerate(["pylimma_logo.png", "omicverse_logo.png"]):
    axl = fig.add_subplot(logo[k])
    axl.imshow(plt.imread(ASSETS / f))
    axl.set_axis_off()
top = row1[1].subgridspec(2, 5, width_ratios=[1, 1, 1, 1, 0.05], wspace=0.12, hspace=0.38)
cats, names = [1, 0, -1], ["up", "ns", "dn"]
for i, (task, coh) in enumerate(zip(full, cohorts)):
    ax = fig.add_subplot(top[i // 4, i % 4])
    cr, cp = calls(load("limma", task), "eb"), calls(load("pylimma", task), "eb")
    M = np.array([[((cr == a) & (cp == b)).sum() for b in cats] for a in cats], float)
    P = M / np.maximum(M.sum(1, keepdims=True), 1) * 100
    im = ax.imshow(P, cmap=SEQ, vmin=0, vmax=100)
    for a in range(3):
        for b in range(3):
            ax.text(b, a, f"{int(M[a, b])}", ha="center", va="center", fontsize=4.5, color="white" if P[a, b] > 55 else INK)
    ax.set_xticks(range(3), names if i // 4 == 1 else [""] * 3, fontsize=5)
    ax.set_yticks(range(3), names if i % 4 == 0 else [""] * 3, fontsize=5)
    ax.tick_params(length=0, pad=1)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f"{coh} (n={int(tab.set_index('task').loc[task, 'samples'])})", fontsize=6, pad=2)
    if i % 4 == 0:
        ax.set_ylabel("limma", fontsize=6, labelpad=1)
    if i // 4 == 1:
        ax.set_xlabel("pylimma", fontsize=6, labelpad=1)
    if i == 0:
        panel(ax, "a", dx=-0.32)
cax = fig.add_subplot(top[:, 4])
cb = fig.colorbar(im, cax=cax)
cb.set_label("% of limma row", fontsize=5.5)
cb.ax.tick_params(labelsize=5)
cb.outline.set_visible(False)

mid = gs[1].subgridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.28)
fg = pd.read_csv(FGSEA, sep="\t")
nesw = fg.pivot_table(index=["task", "pathway"], columns="engine", values="NES")
padj = fg.pivot_table(index=["task", "pathway"], columns="engine", values="padj")
axb = fig.add_subplot(mid[0])
axb.axhline(0, color=INK2, lw=0.6)
axb.scatter(nesw["limma"], nesw["limma_seed2"] - nesw["limma"], s=7, color="#a3a29d", edgecolor="white", linewidth=0.3,
            label="limma, other fgsea seed")
axb.scatter(nesw["limma"], nesw["pylimma"] - nesw["limma"], s=7, color=C_PY, edgecolor="white", linewidth=0.3,
            label="pylimma, same seed")
axb.set_xlabel("NES from limma statistics (Hallmark, 8 cohorts)")
axb.set_ylabel("NES difference from limma")
axb.grid(color=GRID, lw=0.5)
axb.set_axisbelow(True)
dmax_py = (nesw["pylimma"] - nesw["limma"]).abs().max()
dmax_seed = (nesw["limma_seed2"] - nesw["limma"]).abs().max()
axb.legend(fontsize=5, loc="upper left", handletextpad=0.2, borderaxespad=0.3)
axb.text(0.0, -0.75 * dmax_seed, f"max |Δ|\npylimma {dmax_py:.0e}\nother seed {dmax_seed:.2f}",
         ha="center", va="center", fontsize=5, color=INK2)
panel(axb, "b", dx=-0.16)

axc = fig.add_subplot(mid[1])
for k, col in zip(ERRK, QCOL):
    xs = np.sort(np.log10(np.maximum(np.concatenate(errs[k]), 1e-17)))
    axc.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=col, lw=1.2, label=f"{ERRK[k]} (max {10 ** xs[-1]:.0e})")
axc.axvline(-6, color=INK2, lw=0.6, ls=(0, (2, 2)))
axc.set_xlabel("log$_{10}$ |pylimma − limma| / |limma|  (per gene, 8 cohorts)")
axc.set_ylabel("cumulative fraction of genes")
axc.set_xlim(-17, 0)
axc.grid(color=GRID, lw=0.5)
axc.legend(fontsize=4.8, loc="lower right", handlelength=1.4, borderaxespad=0.3)
panel(axc, "c", dx=-0.1)

bot = gs[2].subgridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.24)
axd = fig.add_subplot(bot[0])
x = np.arange(len(cohorts))
for k, (e, col, lab) in enumerate([("limma", C_R, "limma 3.68.5 (R)"), ("pylimma", C_PY, "pylimma")]):
    axd.bar(x + (k - 0.5) * 0.36, tim[e].to_numpy(), 0.33, color=col, label=lab)
axd.set_xticks(x, [f"{c}\n({n})" for c, n in zip(cohorts, tim.samples)], fontsize=5)
axd.set_ylabel("wall time, 1 thread (s)")
axd.grid(axis="y", color=GRID, lw=0.5)
axd.set_axisbelow(True)
axd.legend(fontsize=5, loc="upper right", borderaxespad=0.2)
panel(axd, "d", dx=-0.08)

axe = fig.add_subplot(bot[1])
ss = tab[tab.kind == "sub"]
for k, (t, col) in enumerate(zip(TESTS, [QCOL[0], QCOL[2]])):
    for i, n in enumerate(sorted(ss.n.unique())):
        v = ss.loc[ss.n == n, f"{t}_p_maxrel"].to_numpy()
        jit = np.random.default_rng(i).uniform(-0.07, 0.07, len(v))
        axe.scatter(np.full(len(v), i + (k - 0.5) * 0.3) + jit, np.log10(np.maximum(v, 1e-17)), s=6, color=col,
                    edgecolor="white", linewidth=0.3, label=TESTS[t] if i == 0 else None)
axe.axhline(-6, color=INK2, lw=0.6, ls=(0, (2, 2)))
axe.set_xticks(range(3), [f"{int(n)} vs {int(n)}" for n in sorted(ss.n.unique())])
axe.set_ylabel("log$_{10}$ max rel. difference of P")
axe.set_xlabel("TCGA-LUAD, normal vs tumour (20 draws each)")
axe.grid(axis="y", color=GRID, lw=0.5)
axe.set_ylim(-10.2, -5)
axe.legend(fontsize=5, loc="lower center", ncol=2, columnspacing=0.6, handletextpad=0.1, borderaxespad=0.2)
panel(axe, "e", dx=-0.12)
fig.savefig(OUT / "fig1.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig1.png", bbox_inches="tight", dpi=200)

# ======================= Figure S1 =======================
fig, axs = plt.subplots(2, 4, figsize=(7.2, 3.8))
for ax, task, coh in zip(axs.flat, full, cohorts):
    r, p = load("limma", task), load("pylimma", task)
    a = -np.log10(np.maximum(r.eb_PValue, 1e-320))
    b = -np.log10(np.maximum(p.eb_PValue, 1e-320))
    ax.scatter(a, b, s=1.5, color=C_PY, linewidth=0, rasterized=True)
    lim = max(a.max(), b.max()) * 1.03
    ax.plot([0, lim], [0, lim], color=INK2, lw=0.5)
    ax.set_title(f"TCGA-{coh}", fontsize=6, pad=2)
    ax.text(0.04, 0.9, f"r = {np.corrcoef(a, b)[0, 1]:.12f}", transform=ax.transAxes, fontsize=5)
for ax in axs[1]:
    ax.set_xlabel("limma −log$_{10}$P (eBayes)")
for ax in axs[:, 0]:
    ax.set_ylabel("pylimma −log$_{10}$P")
fig.tight_layout()
fig.savefig(OUT / "figS1_pvalue_scatter.pdf", bbox_inches="tight")
fig.savefig(OUT / "figS1_pvalue_scatter.png", bbox_inches="tight", dpi=200)

# ======================= numbers quoted in the text =======================
F, S = tab[tab.kind == "full"], tab[tab.kind == "sub"]
sig = padj <= 0.05
sp_ = tim.limma / tim.pylimma
macros = {
    "nSamplesTotal": f"{int(F.samples.sum()):,}", "nGenesMin": f"{int(F.genes_kept.min()):,}",
    "nGenesMax": f"{int(F.genes_kept.max()):,}",
    "ebDegTotal": f"{int(F.eb_deg_R.sum()):,}", "trDegTotal": f"{int(F.tr_deg_R.sum()):,}",
    "discordFull": int(F[["eb_discordant", "tr_discordant"]].to_numpy().sum()),
    "fdrFlip": int(F[["eb_fdr_flip", "tr_fdr_flip"]].to_numpy().sum()),
    "discordSub": int(S[["eb_discordant", "tr_discordant"]].to_numpy().sum()),
    "ebPmaxrel": sci(F.eb_p_maxrel.max()), "trPmaxrel": sci(F.tr_p_maxrel.max()),
    "subPmaxrel": sci(S[["eb_p_maxrel", "tr_p_maxrel"]].to_numpy().max()),
    "ebPmedrel": sci(float(np.median(np.concatenate(errs["eb_PValue"])))),
    "wMaxrel": sci(np.concatenate(errs["w_mean"]).max()), "sigmaMaxrel": sci(np.concatenate(errs["sigma"]).max()),
    "lfcMaxrel": sci(np.concatenate(errs["eb_logFC"]).max()),
    "dfpriorMaxrel": sci(tab.dfprior_rel.max()), "s2priorMaxrel": sci(tab.s2prior_rel.max()),
    "fgseaPairs": int(sig["limma"].size), "fgseaSigR": int(sig["limma"].sum()),
    "fgseaDiscordPy": int((sig["limma"] != sig["pylimma"]).sum()),
    "fgseaDiscordSeed": int((sig["limma"] != sig["limma_seed2"]).sum()),
    "fgseaNesPy": sci(dmax_py), "fgseaNesSeed": f"{dmax_seed:.2f}",
    "tBRCAR": f"{tim.loc['BRCA', 'limma']:.0f}", "tBRCAPy": f"{tim.loc['BRCA', 'pylimma']:.0f}",
    "speedMin": f"{sp_.min():.2f}", "speedMax": f"{sp_.max():.2f}",
}
OLD = RUNS / "pylimma_v010"
if OLD.exists():      # the same benchmark with python-limma 0.1.0
    flips, ratios = 0, []
    for task in full:
        r, o = load("limma", task), pd.read_csv(OLD / task / "results.tsv.gz", sep="\t")
        flips += int(((r.eb_FDR <= FDR) != (o.eb_FDR <= FDR)).sum())
        ratios.append(float(np.nanmin(o.eb_t / r.eb_t)))
    macros.update(oldFlip=f"{flips:,}", oldTminRatio=f"{min(ratios):.2f}",
                  oldCohorts=sum(1 for task in full if ((load("limma", task).eb_FDR <= FDR) !=
                                 (pd.read_csv(OLD / task / "results.tsv.gz", sep="\t").eb_FDR <= FDR)).any()))
macros.update(bMaxabs=sci(max(np.nanmax(np.abs(load("limma", t).eb_B - load("pylimma", t).eb_B)) for t in full)))
with open(OUT / "numbers.tex", "w") as fh:
    for k, v in macros.items():
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(open(OUT / "numbers.tex").read())
print(tim.to_string())
print(stages.to_string())

# pylimma technical report: benchmark

The report is published at <https://omicverse.org/articles/pylimma.html>. It compares pylimma with
Bioconductor limma 3.68.5 (limma-voom: `voom` → `lmFit` → `eBayes` / `treat`) on eight TCGA tumour–normal
cohorts and 60 small-sample draws. This directory holds the benchmark that produces its numbers and figures.

## Reproduce

The tasks (TCGA counts and LUAD draws) are built by `py-edgeR/report/scripts/01_prepare_tcga.py`; see
[py-edgeR/report](https://github.com/omicverse/py-edgeR/tree/main/report) for the downloads. Requirements: a
Python env with `python-limma`, `pyedger`, pandas, matplotlib; an R env with limma 3.68.5, edgeR 4.10.5, fgsea,
jsonlite, data.table.

```bash
TASKS=$DATA/edger_tasks OUT=$DATA/limma_runs scripts/run_all.sh          # limma and pylimma, 1 thread each
Rscript scripts/04_fgsea.R $DATA/limma_runs $DATA/edger_tasks/gene_symbols.tsv \
        $DATA/genesets/MSigDB_Hallmark_2020.gmt $DATA/limma_runs/fgsea_hallmark.tsv
python scripts/05_figures.py $DATA/limma_runs $DATA/limma_runs/fgsea_hallmark.tsv $DATA/edger_tasks figures
```

`05_figures.py` writes `figures/numbers.tex`, the numbers quoted in the article.

#!/bin/bash
# Benchmark driver: R limma and pylimma, one process per cohort in parallel, so both are
# timed under the same load. Everything single-threaded (BLAS pinned to one thread).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
TASKS=${TASKS:?set TASKS to the task directory made by py-edgeR/report/scripts/01_prepare_tcga.py}
OUT=${OUT:?set OUT}
RSCRIPT=${RSCRIPT:-Rscript}
PY=${PY:-python}
JOBS=${JOBS:-9}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export HERE TASKS OUT RSCRIPT PY
GLOBS="full_BRCA full_COAD full_KIRC full_LIHC full_LUAD full_LUSC full_PRAD full_THCA sub*"
echo "== limma"; date
printf '%s\n' $GLOBS | xargs -P "$JOBS" -I{} sh -c '"$RSCRIPT" "$HERE/02_run_limma.R" "$TASKS" "$OUT/limma" 1 "{}"'
echo "== pylimma"; date
printf '%s\n' $GLOBS | xargs -P "$JOBS" -I{} sh -c '"$PY" "$HERE/03_run_pylimma.py" "$TASKS" "$OUT/pylimma" 1 "{}"'
date

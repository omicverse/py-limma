# Reference run: Bioconductor limma-voom on every benchmark task.
# Usage: Rscript 02_run_limma.R <tasks_dir> <out_dir> <n_timing_reps> [task_glob]
# Pipeline (limma 3.68 defaults, as in the limma user's guide):
#   edgeR filterByExpr -> normLibSizes (TMM) -> voom(design) -> lmFit
#   -> eBayes (moderated t)   and   -> treat(lfc = 1)
suppressMessages({library(edgeR); library(limma); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
tasks_dir <- args[1]; out_dir <- args[2]; reps <- as.integer(args[3])
pattern <- if (length(args) >= 4) args[4] else "*"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

run_pipeline <- function(X, group) {
  tm <- list()
  t0 <- proc.time()[["elapsed"]]
  y <- DGEList(X, group = group)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  design <- model.matrix(~group)
  tm$filter_norm <- proc.time()[["elapsed"]] - t0; t0 <- proc.time()[["elapsed"]]
  v <- voom(y, design)
  tm$voom <- proc.time()[["elapsed"]] - t0; t0 <- proc.time()[["elapsed"]]
  fit <- lmFit(v, design)
  efit <- eBayes(fit)
  tm$lmFit_eBayes <- proc.time()[["elapsed"]] - t0; t0 <- proc.time()[["elapsed"]]
  tfit <- treat(fit, lfc = 1)
  tm$treat <- proc.time()[["elapsed"]] - t0
  list(v = v, keep = keep, efit = efit, tfit = tfit, time = tm)
}

tasks <- Sys.glob(file.path(tasks_dir, pattern))
tasks <- tasks[file.exists(file.path(tasks, "counts.tsv.gz"))]
for (td in tasks) {
  name <- basename(td)
  od <- file.path(out_dir, name)
  if (file.exists(file.path(od, "timing.json"))) next
  dir.create(od, showWarnings = FALSE)
  X <- as.matrix(read.delim(gzfile(file.path(td, "counts.tsv.gz")), row.names = 1, check.names = FALSE))
  s <- read.delim(file.path(td, "samples.tsv"), row.names = 1)
  group <- factor(s[colnames(X), "group"], levels = c("Normal", "Tumor"))
  times <- list()
  for (r in seq_len(reps)) {
    res <- run_pipeline(X, group)
    times[[r]] <- res$time
  }
  e <- res$efit; t <- res$tfit; v <- res$v
  out <- data.frame(
    gene = rownames(v),
    AveExpr = e$Amean,
    w_mean = rowMeans(v$weights),
    w_first = v$weights[, 1],
    E_first = v$E[, 1],
    sigma = e$sigma,
    s2post = e$s2.post,
    eb_logFC = e$coefficients[, 2],
    eb_t = e$t[, 2],
    eb_PValue = e$p.value[, 2],
    eb_FDR = p.adjust(e$p.value[, 2], "BH"),
    eb_B = e$lods[, 2],
    tr_t = t$t[, 2],
    tr_PValue = t$p.value[, 2],
    tr_FDR = p.adjust(t$p.value[, 2], "BH")
  )
  write.table(out, gzfile(file.path(od, "results.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
  scal <- list(n_keep = sum(res$keep), norm_factors = unname(v$targets$norm.factors),
               df_prior = e$df.prior, s2_prior = e$s2.prior)
  write_json(scal, file.path(od, "scalars.json"), digits = NA, auto_unbox = TRUE)
  write_json(list(engine = "limma", version = as.character(packageVersion("limma")),
                  R = R.version.string, reps = times), file.path(od, "timing.json"), digits = NA, auto_unbox = TRUE)
  cat(name, sprintf("%.2fs", sum(unlist(times[[reps]]))), "\n")
}

# Pre-ranked GSEA (fgsea, MSigDB Hallmark 2020) on the moderated t of every full cohort,
# once from limma's and once from pylimma's statistics, with the same seed.
# A second limma run with another seed measures fgsea's own Monte-Carlo noise.
# Usage: Rscript 04_fgsea.R <runs_dir> <symbols.tsv> <hallmark.gmt> <out.tsv>
suppressMessages({library(fgsea); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
runs <- args[1]; sym <- fread(args[2]); gmt <- gmtPathways(args[3]); out <- args[4]
sym[, ens := sub("\\..*$", "", id)]

rank_stat <- function(f) {
  r <- as.data.table(read.delim(gzfile(f)))
  r[, ens := sub("\\..*$", "", gene)]
  r <- merge(r, sym[, .(ens, gene_symbol = gene)], by = "ens")
  r[, stat := eb_t]
  r <- r[order(-abs(stat))][!duplicated(gene_symbol)]
  setNames(r$stat, r$gene_symbol)
}

res <- list()
for (td in Sys.glob(file.path(runs, "limma", "full_*"))) {
  task <- basename(td)
  for (run in list(c("limma", "limma", 1), c("pylimma", "pylimma", 1), c("limma_seed2", "limma", 2))) {
    st <- rank_stat(file.path(runs, run[2], task, "results.tsv.gz"))
    set.seed(as.integer(run[3]))
    fg <- fgsea(gmt, st, minSize = 15, maxSize = 500, eps = 0, nproc = 1)
    res[[length(res) + 1]] <- data.table(task = task, engine = run[1], pathway = fg$pathway,
                                         NES = fg$NES, pval = fg$pval, padj = fg$padj, size = fg$size)
  }
  cat(task, "\n")
}
fwrite(rbindlist(res), out, sep = "\t")

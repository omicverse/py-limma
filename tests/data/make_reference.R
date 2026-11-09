# Regenerate the limma reference values used by tests/test_r_parity.py.
# Run from the repository root:  Rscript tests/data/make_reference.R
# Reference: R 4.5.3, limma 3.68.5, edgeR 4.10.5.
# counts.tsv: 600 simulated genes, 6 libraries x {NT, KO} (the py-edgeR test counts).
suppressMessages({library(edgeR); library(limma)})
d <- "tests/data"
w <- function(x, f) write.table(x, file.path(d, f), sep = "\t", quote = FALSE, row.names = FALSE)
g17 <- function(x) sprintf("%.17g", x)
X <- as.matrix(read.delim(file.path(d, "counts.tsv"), row.names = 1))
lib <- factor(sub("_(KO|NT)$", "", colnames(X)))
cond <- factor(sub("^.*_", "", colnames(X)), levels = c("NT", "KO"))
y <- DGEList(X, group = cond); keep <- filterByExpr(y); y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
w(data.frame(keep = which(keep) - 1L, lib_size = NA)[, 1, drop = FALSE], "keep.tsv")
w(data.frame(lib = g17(y$samples$lib.size * y$samples$norm.factors)), "eff_lib.tsv")

# 1) voom + lmFit + eBayes, design ~ lib + cond
D1 <- model.matrix(~ lib + cond)
v <- voom(y, D1)
w(data.frame(E = g17(as.vector(v$E)), W = g17(as.vector(v$weights))), "voom.tsv")
f1 <- eBayes(lmFit(v, D1))
k <- ncol(D1)
w(data.frame(coef = g17(f1$coefficients[, k]), t = g17(f1$t[, k]), p = g17(f1$p.value[, k]), B = g17(f1$lods[, k]),
             s2post = g17(f1$s2.post)), "ebayes.tsv")
w(data.frame(df_prior = g17(f1$df.prior), s2_prior = g17(f1$s2.prior), var_prior = g17(f1$var.prior[k])), "ebayes_scalar.tsv")
tt <- topTable(f1, coef = k, number = Inf)
w(data.frame(gene = rownames(tt), B = g17(tt$B), adj = g17(tt$adj.P.Val)), "toptable.tsv")

# 2) a continuous covariate: its stdev.unscaled is far below 0.1, so clipping it (pyedger-limma <= 0.1.0) deflated t
age <- seq(25, 80, length.out = ncol(X))
D2 <- cbind(Intercept = 1, KO = as.numeric(cond == "KO"), age = age)
f2 <- eBayes(lmFit(voom(y, D2), D2))
w(data.frame(t = g17(f2$t[, 3]), p = g17(f2$p.value[, 3]), B = g17(f2$lods[, 3]), su = g17(f2$stdev.unscaled[, 3])), "ebayes_age.tsv")
write.table(data.frame(age = g17(age)), file.path(d, "age.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 3) treat, lfc = 1
tr <- treat(lmFit(v, D1), lfc = 1)
w(data.frame(t = g17(tr$t[, k]), p = g17(tr$p.value[, k])), "treat.tsv")

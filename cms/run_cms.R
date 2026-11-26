#!/usr/bin/env Rscript
# CMS calls (CMSclassifier random forest) vs mutation clusters. Run from anywhere:
#   Rscript cms/run_cms.R          full analysis (needs CMSclassifier in cms/rlib)
#   Rscript cms/run_cms.R --prep   only build and check the input
# The classifier wants log2 expression, Entrez rows, and centres each gene over the cohort,
# so gene length cancels and log2(CPM + 1) is used. Centring uses all 577 selected tumours;
# proportions are reported for the 552 in C1/C2.
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
root <- dirname(here)
.libPaths(c(file.path(here, "rlib"), .libPaths()))
suppressPackageStartupMessages(library(org.Hs.eg.db))
out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE)

assign <- read.delim(file.path(root, "results/02_clustering/cluster_assignments.tsv"))
meta <- read.delim(file.path(root, "results/01_mutation_matrix/sample_annotation.tsv"))
read_expr <- function(p) {
  x <- as.matrix(read.delim(gzfile(p), row.names = 1, check.names = FALSE))
  x[!grepl("_PAR_Y$", rownames(x)), intersect(colnames(x), assign$sample), drop = FALSE]
}
raw <- file.path(root, "data/raw")
a <- read_expr(file.path(raw, "TCGA-COAD.star_counts.tsv.gz"))
b <- read_expr(file.path(raw, "TCGA-READ.star_counts.tsv.gz"))
stopifnot(setequal(rownames(a), rownames(b)))
counts <- pmax(round(2^cbind(a, b[rownames(a), ]) - 1), 0)[, assign$sample]
lib <- colSums(counts)
ens <- sub("\\..*", "", rownames(counts))
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, unique(ens), "ENTREZID", "ENSEMBL"))
map <- map[!is.na(map$ENTREZID), ]
map <- map[!map$ENSEMBL %in% map$ENSEMBL[duplicated(map$ENSEMBL)], ]  # ambiguous IDs
keep <- ens %in% map$ENSEMBL
entrez <- rowsum(counts[keep, ], map$ENTREZID[match(ens[keep], map$ENSEMBL)])
expr <- as.data.frame(log2(t(t(entrez) / lib) * 1e6 + 1), check.names = FALSE)
cat("input:", nrow(expr), "Entrez genes x", ncol(expr), "tumours\n")
if ("--prep" %in% args) quit(save = "no")

suppressPackageStartupMessages(library(CMSclassifier))
genes <- listModelGenes("RF")
cat("model genes:", length(genes), "missing from input:", length(setdiff(genes, rownames(expr))), "\n")
rf <- classifyCMS.RF(expr)
write.table(cbind(sample = rownames(rf), rf), file.path(out, "cms_rf_calls.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

lv <- c("CMS1", "CMS2", "CMS3", "CMS4", "Unclassified")
d <- data.frame(assign, cms = rf[assign$sample, "RF.predictedCMS"],
                hypermutated = meta$hypermutated[match(assign$sample, meta$sample)])
d$cms[is.na(d$cms) | !d$cms %in% lv] <- "Unclassified"   # posterior < 0.5
d$cms <- factor(d$cms, lv)
write.table(d, file.path(out, "cms_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
d <- d[d$cluster %in% c("C1", "C2"), ]; d$cluster <- factor(d$cluster)
tab <- table(d$cluster, d$cms); prop <- prop.table(tab, 1)
long <- as.data.frame(tab); names(long) <- c("cluster", "cms", "n")
long$proportion <- as.data.frame(prop)$Freq
write.table(long, file.path(out, "cms_proportions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tab); print(round(100 * prop, 1))
set.seed(501)
pv <- function(t) fisher.test(t, simulate.p.value = TRUE, B = 1e5)$p.value
called <- tab[, lv[1:4]]
cat("Fisher (simulated) p, classified only:", pv(called[, colSums(called) > 0]),
    " incl. unclassified:", pv(tab[, colSums(tab) > 0]), "\n")
print(table(d$cluster, d$cms, d$hypermutated))

# order and colours: CMS convention (Guinney 2015); grey = unclassified
cols <- c("#E69F00", "#0072B2", "#CC79A7", "#009E73", "#BDBDBD")
png(file.path(out, "cms_proportions.png"), 1500, 1100, res = 200)
par(mar = c(4, 4.5, 2.5, 8), xpd = NA)
bp <- barplot(100 * t(prop), col = cols, border = "white", las = 1, space = 0.5,
  ylab = "Tumours (%)", main = "CMS composition of mutation clusters",
  names.arg = sprintf("%s (n=%d)", rownames(tab), rowSums(tab)))
for (i in 1:2) { p <- 100 * prop[i, ]; mid <- cumsum(p) - p / 2; s <- p >= 4
  if (any(s)) text(bp[i], mid[s], sprintf("%.0f%%", p[s])) }
legend(max(bp) + 0.7, 100, rev(lv), fill = rev(cols), border = NA, bty = "n")
invisible(dev.off())

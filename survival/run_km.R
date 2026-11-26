#!/usr/bin/env Rscript
# Overall survival, C1 vs C2: Kaplan-Meier, log-rank, Cox.   Rscript survival/run_km.R
# Time = days to death (dead) or days to last follow-up (alive); time <= 0 or missing is excluded.
library(survival)
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
root <- dirname(here); out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE)
cl <- read.delim(file.path(root, "results/02_clustering/cluster_assignments.tsv"))
ann <- read.delim(file.path(root, "results/01_mutation_matrix/sample_annotation.tsv"))
want <- c("sample", "vital_status.demographic", "days_to_death.demographic",
          "days_to_last_follow_up.diagnoses", "age_at_index.demographic")
clin <- do.call(rbind, lapply(c("COAD", "READ"), function(p) read.delim(gzfile(
  file.path(root, "data/raw", sprintf("TCGA-%s.clinical.tsv.gz", p))), quote = "", check.names = FALSE)[, want]))
names(clin) <- c("sample", "vital", "death", "fu", "age")
clin <- clin[!duplicated(clin$sample), ]
d <- merge(merge(cl[cl$cluster %in% c("C1", "C2"), ], clin, all.x = TRUE), ann[, c("sample", "stage")])
d$event <- ifelse(d$vital == "Dead", 1, ifelse(d$vital == "Alive", 0, NA))
d$days <- ifelse(d$event == 1, d$death, d$fu)
d$years <- d$days / 365.25; d$cluster <- factor(d$cluster)
cat("clustered", nrow(d), "| missing time/status", sum(is.na(d$days)), "| time <= 0:",
    sum(d$days <= 0, na.rm = TRUE), "of which deaths", sum(d$days <= 0 & d$event == 1, na.rm = TRUE), "\n")
k <- d[!is.na(d$days) & d$days > 0, ]
write.table(k[, c("sample", "cluster", "days", "event", "age", "stage")],
            file.path(out, "km_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
fit <- survfit(Surv(years, event) ~ cluster, k); print(fit)
print(summary(fit, times = c(1, 3, 5), extend = TRUE))
lr <- survdiff(Surv(years, event) ~ cluster, k); p <- pchisq(lr$chisq, 1, lower.tail = FALSE); print(lr)
cx <- coxph(Surv(years, event) ~ cluster, k); print(summary(cx)); print(cox.zph(cx))
k$stage <- factor(k$stage); print(table(k$stage, k$cluster, useNA = "ifany"))
s <- droplevels(k[k$stage %in% c("I", "II", "III", "IV") & !is.na(k$age), ])
cat("adjusted for age + stage, n =", nrow(s), "\n")
ca <- coxph(Surv(years, event) ~ cluster + age + stage, s); print(summary(ca)); print(cox.zph(ca))
z <- d[!is.na(d$days), ]; z$years <- pmax(z$years, 0.5 / 365.25)
cat("SENSITIVITY incl. zero times, n =", nrow(z), "\n")
print(survdiff(Surv(years, event) ~ cluster, z)); print(summary(coxph(Surv(years, event) ~ cluster, z))$conf.int)

cols <- c("#6250d6", "#e87ba4"); tk <- 0:8
png(file.path(out, "km_clusters.png"), 1500, 1300, res = 200)
par(mar = c(9, 7, 2.5, 1.5), las = 1, xpd = NA)
plot(NA, xlim = c(0, 8), ylim = c(0, 1), xaxs = "i", yaxs = "i", xaxt = "n", bty = "l",
     xlab = "", ylab = "Overall survival probability", main = "Overall survival by mutation cluster")
axis(1, tk); mtext("Years", 1, 2.2)
par(xpd = FALSE)
for (i in 1:2) { f <- fit[i]
  x <- c(0, f$time); r <- function(v) c(v[1], rep(v, each = 2)[-c(1)])
  ok <- !is.na(f$lower) & !is.na(f$upper)
  xs <- c(0, rep(f$time[ok], each = 2)); lo <- c(1, 1, rep(f$lower[ok], each = 2)); hi <- c(1, 1, rep(f$upper[ok], each = 2))
  m <- length(xs); polygon(c(xs, rev(xs)), c(lo[1:m], rev(hi[1:m])), col = adjustcolor(cols[i], .15), border = NA) }
lines(fit, col = cols, lwd = 2, lty = c(1, 2), mark.time = TRUE, conf.int = FALSE, cex = .5)
ci <- summary(cx)$conf.int
legend("bottomleft", sprintf("%s (n=%d, %d deaths)", levels(k$cluster), table(k$cluster), tapply(k$event, k$cluster, sum)),
       col = cols, lty = 1:2, lwd = 2, bty = "n", cex = .8, inset = c(.02, .14))
legend("bottomleft", c(sprintf("Log-rank p = %.2g", p), sprintf("HR C2 vs C1 = %.2f (95%% CI %.2f-%.2f)", ci[1], ci[3], ci[4])),
       bty = "n", cex = .8, inset = c(.02, .01))
risk <- sapply(tk, function(t) tapply(k$years >= t, k$cluster, sum))
mtext("Number at risk", 1, 4, at = 0, adj = 0, cex = .8, font = 2)
for (i in 1:2) { mtext(levels(k$cluster)[i], 1, 4 + i * 1.1, at = -0.6, cex = .8, col = "black", font = 2)
  mtext(risk[i, ], 1, 4 + i * 1.1, at = tk, cex = .8) }
invisible(dev.off())

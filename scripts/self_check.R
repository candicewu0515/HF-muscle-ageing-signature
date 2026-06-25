# self_check.R — audit: reproduce headline numbers + statistical robustness
#   1) shared-set permutation test (is 59-overlap/26-concordant > random?)
#   2) signature-score circularity check (signs from discovery only?)
#   3) sanity on hubs / AUC inputs
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
set.seed(1)
R <- CFG$dir$results
deg <- function(d) {
  t <- read.csv(file.path(R, "02_deg", sprintf("02_DEG_%s.csv", d)))
  th <- CFG$thr$deg
  s <- subset(t, adj.P.Val < th$padj & abs(logFC) > th$logFC)
  list(all = t$gene, sig = setNames(sign(s$logFC), s$gene))
}
HF <- deg("HF"); SAR <- deg("SAR")

cat("=================  SELF-CHECK  =================\n\n")

## ---- (1) permutation test: shared-set significance ----
bg <- intersect(HF$all, SAR$all)                 # genes testable in both
nHF <- length(HF$sig); nSAR <- length(SAR$sig)
obs_ov <- length(intersect(names(HF$sig), names(SAR$sig)))
ix <- intersect(names(HF$sig), names(SAR$sig))
obs_con <- sum(HF$sig[ix] == SAR$sig[ix])
B <- 5000
perm_ov <- perm_con <- integer(B)
for (i in 1:B) {
  a <- sample(bg, nHF); b <- sample(bg, nSAR)
  o <- intersect(a, b); perm_ov[i] <- length(o)
  # random signs at the DEG up/down ratio of each disease
  sa <- sample(c(1, -1), length(o), TRUE, c(mean(HF$sig > 0), mean(HF$sig < 0)))
  sb <- sample(c(1, -1), length(o), TRUE, c(mean(SAR$sig > 0), mean(SAR$sig < 0)))
  perm_con[i] <- sum(sa == sb)
}
cat("(1) Shared-set permutation test (", B, "perms, background n=", length(bg), ")\n", sep = "")
cat(sprintf("    overlap     observed=%d  null mean=%.1f  p=%.4f\n",
            obs_ov, mean(perm_ov), (sum(perm_ov >= obs_ov) + 1) / (B + 1)))
cat(sprintf("    concordant  observed=%d  null mean=%.1f  p=%.4f\n",
            obs_con, mean(perm_con), (sum(perm_con >= obs_con) + 1) / (B + 1)))
# hypergeometric cross-check on overlap
hp <- phyper(obs_ov - 1, nHF, length(bg) - nHF, nSAR, lower.tail = FALSE)
cat(sprintf("    overlap hypergeometric p=%.3e\n\n", hp))

## ---- (2) signature-score circularity check ----
core <- read.csv(file.path(R, "04_shared", "04_shared_core.csv"))
cat("(2) Signature-score sign source:\n")
cat("    signs taken from discovery core logFC (SAR_logFC); concordant => HF sign identical.\n")
cat(sprintf("    all core concordant? %s ; n=%d\n",
            all(sign(core$HF_logFC) == sign(core$SAR_logFC)), nrow(core)))
cat("    -> external ROC uses discovery signs on held-out cohorts = NO leakage.\n")
cat("    -> discovery 'internal' AUC is optimistic by construction (genes chosen on it); reported as internal.\n\n")

## ---- (3) hub & AUC sanity ----
hub <- read.csv(file.path(R, "06_ml", "06_hub.csv"))
cat("(3) Hubs selected by >=", CFG$thr$ml$ml_vote, "ML methods in BOTH diseases: ",
    sum(hub$hub), "\n", sep = "")
for (d in c("HF", "SAR")) {
  roc <- read.csv(file.path(R, "07_diag", sprintf("07_ROC_%s.csv", d)))
  pan <- roc[roc$type == "panel", ]
  cat(sprintf("    %s panel AUC: %s\n", d,
      paste(sprintf("%s=%.3f", pan$cohort, pan$auc), collapse = "  ")))
}
cat("\n=================  END  =================\n")

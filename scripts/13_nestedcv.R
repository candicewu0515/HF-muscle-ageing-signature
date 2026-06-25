# =====================================================================
# 13_nestedcv.R — honest out-of-sample AUC for the signature score
#   Addresses circularity: the directional signature weights (sign of logFC)
#   are RE-DERIVED inside each training fold; the held-out fold is scored
#   with those signs and z-scored within itself (unsupervised) -> no leakage.
#   Replaces the optimistic "internal AUC" (1.00) with cross-validated AUC.
#   output: results/13_robust/13_nestedcv_auc.csv
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(pROC); library(limma) })
set.seed(CFG$SEED)
qc <- file.path(CFG$dir$results, "01_qc")
out <- file.path(CFG$dir$results, "13_robust"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
hub <- read.csv(file.path(CFG$dir$results, "06_ml", "06_hub.csv")); hub <- hub$gene[hub$hub]

fold_signs <- function(ex, y) {                 # sign of case-vs-control logFC per gene
  d <- model.matrix(~ y); fit <- eBayes(lmFit(ex, d))
  sign(topTable(fit, coef = 2, number = Inf, sort.by = "none")$logFC)
}

nested_auc <- function(dis, K = 5, reps = 20) {
  acc <- CFG$datasets[[dis]]$discovery$acc
  d <- readRDS(file.path(qc, sprintf("01_expr_%s_%s.rds", dis, acc)))
  g <- intersect(hub, rownames(d$expr)); ex <- d$expr[g, ]; y <- d$group
  aucs <- c()
  for (r in seq_len(reps)) {
    set.seed(100 + r)
    folds <- sample(rep(seq_len(K), length.out = ncol(ex)))
    oof_score <- rep(NA_real_, ncol(ex))
    for (k in seq_len(K)) {
      tr <- folds != k; te <- !tr
      if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) next
      sg <- fold_signs(ex[, tr, drop = FALSE], droplevels(y[tr]))     # signs from TRAIN only
      Zte <- t(scale(t(ex[, te, drop = FALSE])))                       # z within test (unsupervised)
      Zte[is.na(Zte)] <- 0
      oof_score[te] <- as.numeric(t(Zte) %*% sg)
    }
    ok <- !is.na(oof_score)
    aucs <- c(aucs, as.numeric(roc(y[ok], oof_score[ok], quiet = TRUE, direction = "<")$auc))
  }
  data.frame(disease = dis, cohort = paste0(acc, " (nested ", K, "-fold x", reps, ")"),
             auc_mean = mean(aucs), auc_sd = sd(aucs),
             auc_lo = quantile(aucs, .025), auc_hi = quantile(aucs, .975),
             n = ncol(ex), genes = length(g))
}

res <- do.call(rbind, lapply(active_diseases(), nested_auc))
write.csv(res, file.path(out, "13_nestedcv_auc.csv"), row.names = FALSE)
for (i in seq_len(nrow(res))) log_msg(res$disease[i], " nested-CV AUC=",
  signif(res$auc_mean[i], 3), " (", signif(res$auc_lo[i], 3), "-", signif(res$auc_hi[i], 3),
  ") vs optimistic internal", step = "13")
log_msg("13 nested-CV done.", step = "13")

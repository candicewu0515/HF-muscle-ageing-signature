# =====================================================================
# 15_meta_auc.R — random-effects meta-analysis of the muscle signature AUC
#   across 5 ageing-muscle cohorts -> pooled AUC + forest plot.
#   Strengthens the (weak, single-cohort) sarcopenia validation.
#   output: results/15_meta/15_muscle_meta_auc.csv  forest_muscle_auc.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(GEOquery); library(limma); library(pROC); library(metafor) })
options(timeout = 1800); Sys.setenv(VROOM_CONNECTION_SIZE = 5e6)
set.seed(CFG$SEED)
R <- CFG$dir$results
out <- file.path(R, "15_meta"); dir.create(out, showWarnings = FALSE, recursive = TRUE)

core <- read.csv(file.path(R, "04_shared", "04_shared_core.csv"))
SGN  <- setNames(sign(core$SAR_logFC), core$gene)
hub  <- read.csv(file.path(R, "06_ml", "06_hub.csv")); hub <- hub$gene[hub$hub]

# compact GPL570 loader (Gene Symbol col); returns gene x sample matrix
load_gpl570 <- function(acc) {
  es <- getGEO(acc, destdir = CFG$dir$geo, GSEMatrix = TRUE, getGPL = TRUE)
  es <- if (is.list(es)) es[[1]] else es
  ex <- Biobase::exprs(es)
  qx <- as.numeric(quantile(ex, c(0, .5, .99, 1), na.rm = TRUE))
  if (isTRUE(qx[3] > 100)) { ex[ex < 0] <- NA; ex <- log2(ex + 1) }
  fd <- Biobase::fData(es)
  sym <- fd[[intersect(c("Gene Symbol", "Gene symbol"), colnames(fd))[1]]]
  sym <- trimws(sub(" ?///.*$", "", as.character(sym)))
  keep <- !is.na(sym) & sym != "" & rowSums(!is.na(ex)) > 0
  ex <- ex[keep, ]; sym <- sym[keep]
  o <- order(rowMeans(ex, na.rm = TRUE), decreasing = TRUE)
  ex <- ex[o, ]; sym <- sym[o]; ex <- ex[!duplicated(sym), ]; rownames(ex) <- sym[!duplicated(sym)]
  list(ex = normalizeBetweenArrays(ex, method = "quantile"), pd = Biobase::pData(es))
}

# old=case / young=control grouping per cohort
group_muscle <- function(acc, pd) {
  g <- switch(acc,
    "GSE9103" = { v <- tolower(pd[["source_name_ch1"]])
      ifelse(grepl("sedentary", v) & grepl("old", v), "case",
        ifelse(grepl("sedentary", v) & grepl("young", v), "control", NA)) },  # sedentary only
    "GSE38718" = { v <- tolower(pd[["age group:ch1"]])
      ifelse(grepl("old", v), "case", ifelse(grepl("young", v), "control", NA)) },
    "GSE28392" = { v <- tolower(pd[["source_name_ch1"]])
      ifelse(grepl("pre", v) & grepl("old", v), "case",
        ifelse(grepl("pre", v) & grepl("young", v), "control", NA)) })           # baseline only
  factor(g, levels = c("control", "case"))
}

sig_score <- function(ex, genes) {
  g <- intersect(genes, intersect(rownames(ex), names(SGN)))
  Z <- t(scale(t(ex[g, , drop = FALSE]))); Z[is.na(Z)] <- 0
  as.numeric(t(Z) %*% SGN[g])
}
auc_ci <- function(y, s) { r <- roc(y, s, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r, method = "delong")); c(auc = as.numeric(r$auc), lo = ci[1], hi = ci[3]) }

rows <- list()
# existing processed cohorts
for (acc in c("GSE1428", "GSE25941")) {
  d <- readRDS(file.path(R, "01_qc", sprintf("01_expr_SAR_%s.rds", acc)))
  y <- as.integer(d$group) - 1; a <- auc_ci(y, sig_score(d$expr, hub))
  rows[[acc]] <- data.frame(cohort = acc, n = length(y), auc = a["auc"], lo = a["lo"], hi = a["hi"])
  log_msg(acc, " AUC=", signif(a["auc"], 3), " n=", length(y), step = "15")
}
# new cohorts
for (acc in c("GSE9103", "GSE38718", "GSE28392")) {
  r <- tryCatch({
    L <- load_gpl570(acc); grp <- group_muscle(acc, L$pd)
    keep <- !is.na(grp); ex <- L$ex[, keep]; y <- as.integer(droplevels(grp[keep])) - 1
    a <- auc_ci(y, sig_score(ex, hub))
    log_msg(acc, " AUC=", signif(a["auc"], 3), " n=", length(y),
            " (case=", sum(y == 1), " ctrl=", sum(y == 0), ")", step = "15")
    data.frame(cohort = acc, n = length(y), auc = a["auc"], lo = a["lo"], hi = a["hi"])
  }, error = function(e) { log_msg(acc, " FAILED: ", conditionMessage(e), step = "15"); NULL })
  if (!is.null(r)) rows[[acc]] <- r
}

df <- do.call(rbind, rows)
write.csv(df, file.path(out, "15_muscle_meta_auc.csv"), row.names = FALSE)

# random-effects meta on logit(AUC); SE from DeLong CI
df$yi <- qlogis(pmin(pmax(df$auc, 1e-3), 1 - 1e-3))
df$sei <- (qlogis(pmin(df$hi, .999)) - qlogis(pmax(df$lo, .001))) / (2 * 1.96)
m <- rma(yi = yi, sei = sei, data = df, method = "REML")
pooled <- plogis(c(m$beta, m$ci.lb, m$ci.ub))
log_msg("MUSCLE meta-AUC (random-effects, ", nrow(df), " cohorts, N=", sum(df$n),
        "): pooled=", signif(pooled[1], 3), " (", signif(pooled[2], 3), "-", signif(pooled[3], 3),
        ") I2=", round(m$I2), "%", step = "15")

# forest shows the per-cohort DeLong CIs (AUC scale) directly from the table, so
# displayed intervals match 15_muscle_meta_auc.csv exactly; pooled = meta diamond.
png(file.path(out, "forest_muscle_auc.png"), 1150, 720, res = 130)
forest(x = df$auc, ci.lb = df$lo, ci.ub = df$hi, refline = 0.5,
       slab = paste0(df$cohort, " (n=", df$n, ")"),
       xlab = "AUC (directional signature score, ageing muscle)",
       xlim = c(-0.6, 1.9), alim = c(0.3, 1.0), at = c(.3,.5,.7,.9,1.0),
       header = c("Cohort", "AUC [95% CI, DeLong]"), col = "steelblue", psize = 1.2)
addpoly(pooled[1], ci.lb = pooled[2], ci.ub = pooled[3], row = -1,
        mlab = sprintf("Random-effects pooled (I2=%.0f%%)", m$I2), col = "tomato")
dev.off()
log_msg("15 meta-AUC done.", step = "15")

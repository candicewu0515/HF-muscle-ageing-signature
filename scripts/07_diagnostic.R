# =====================================================================
# 07_diagnostic.R — diagnostic performance of the shared hub panel
#   train logistic (rms::lrm) on discovery; validate on EXTERNAL cohort(s)
#   ROC + 95% CI (per-gene & panel), nomogram, calibration, decision curve
#   outputs per disease:
#     results/07_diag/07_ROC_<DIS>.csv          (AUC + 95% CI)
#     roc_<DIS>_<cohort>.png  nomogram_<DIS>.png  calibration_<DIS>.png  dca_<DIS>.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(pROC); library(rms); library(ggplot2) })
set.seed(CFG$SEED)
qc   <- file.path(CFG$dir$results, "01_qc")
mldir<- file.path(CFG$dir$results, "06_ml")
outdir <- file.path(CFG$dir$results, "07_diag")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

hub <- read.csv(file.path(mldir, "06_hub.csv"))
hub <- hub$gene[hub$hub]
core <- read.csv(file.path(CFG$dir$results, "04_shared", "04_shared_core.csv"))
SGN  <- setNames(sign(core$SAR_logFC), core$gene)   # concordant -> HF sign identical
log_msg("diagnostic panel = ", length(hub), " hub genes", step = "07")

load_cohort <- function(dis, acc) {
  d <- readRDS(file.path(qc, sprintf("01_expr_%s_%s.rds", dis, acc)))
  g <- intersect(hub, rownames(d$expr))
  X <- as.data.frame(t(d$expr[g, , drop = FALSE]))
  X$.y <- as.integer(factor(d$group, levels = c("control", "case"))) - 1
  list(X = X, genes = g)
}
# directional signature score: sum of z-scored expression * discovery direction
# (transfers across platforms better than a transferred multivariate logistic)
sig_score <- function(X, genes) {
  g <- intersect(genes, names(SGN))
  Z <- scale(X[, g, drop = FALSE])
  as.numeric(Z %*% SGN[g])
}

# manual decision-curve net benefit
dca_df <- function(y, p, model = "panel") {
  th <- seq(0.01, 0.99, by = 0.01); N <- length(y); prev <- mean(y)
  nb <- sapply(th, function(t) {
    tp <- sum(p >= t & y == 1); fp <- sum(p >= t & y == 0)
    (tp / N) - (fp / N) * (t / (1 - t))
  })
  nb_all <- sapply(th, function(t) prev - (1 - prev) * (t / (1 - t)))
  rbind(
    data.frame(thr = th, nb = nb, strat = model),
    data.frame(thr = th, nb = nb_all, strat = "treat all"),
    data.frame(thr = th, nb = 0, strat = "treat none"))
}

run_dx <- function(dis) {
  disc_acc <- CFG$datasets[[dis]]$discovery$acc
  vaccs <- vapply(CFG$datasets[[dis]]$validation, function(v) v$acc, "")
  tr <- load_cohort(dis, disc_acc)
  # panel = hubs present in discovery AND every validation cohort (one model everywhere)
  panel <- tr$genes
  for (acc in vaccs) panel <- intersect(panel, load_cohort(dis, acc)$genes)
  if (length(panel) < 2) { log_msg(dis, " <2 common hubs across cohorts, skip", step = "07"); return() }
  tr$X <- tr$X[, c(panel, ".y")]
  log_msg(dis, " common panel (", length(panel), "/", length(hub), " hubs): ",
          paste(panel, collapse = ","), step = "07")
  assign("dd", datadist(tr$X), envir = .GlobalEnv); options(datadist = "dd")
  form <- as.formula(paste(".y ~", paste(sprintf("`%s`", panel), collapse = " + ")))
  fit <- lrm(form, data = tr$X, x = TRUE, y = TRUE)

  # ---- nomogram ----
  png(file.path(outdir, sprintf("nomogram_%s.png", dis)), 1500, 950, res = 150)
  plot(nomogram(fit, fun = plogis, funlabel = sprintf("P(%s case)", dis)))
  dev.off()
  # ---- calibration (bootstrap on discovery) ----
  cal <- tryCatch(calibrate(fit, method = "boot", B = 200), error = function(e) NULL)
  if (!is.null(cal)) {
    png(file.path(outdir, sprintf("calibration_%s.png", dis)), 950, 950, res = 150)
    plot(cal, main = sprintf("%s calibration (discovery, 200x boot)", dis)); dev.off()
  }

  # ---- signature-score probability model on discovery (for DCA/calibration of score) ----
  tr_score <- sig_score(tr$X, panel)
  sfit <- glm(tr$X$.y ~ tr_score, family = binomial)

  # ---- external validation ROC (directional signature score) ----
  # internal (discovery) ROC first
  rocrows <- list()
  ri <- roc(tr$X$.y, tr_score, quiet = TRUE, direction = "<")
  cii <- as.numeric(ci.auc(ri, conf.level = CFG$thr$roc$ci, method = "delong"))
  rocrows[["discovery"]] <- data.frame(disease = dis, cohort = paste0(disc_acc, "(internal)"),
    type = "panel", n = nrow(tr$X), auc = as.numeric(ri$auc),
    ci_low = cii[1], ci_high = cii[3], genes = length(panel))

  for (acc in vaccs) {
    va <- load_cohort(dis, acc)
    gp <- intersect(panel, va$genes)
    if (length(gp) < 2) { log_msg(dis, "/", acc, " <2 hubs, skip ROC", step = "07"); next }
    sc <- sig_score(va$X, gp)
    lp <- as.numeric(predict(sfit, newdata = data.frame(tr_score = sc), type = "response"))
    y <- va$X$.y
    r <- roc(y, sc, quiet = TRUE, direction = "<")
    ci <- as.numeric(ci.auc(r, conf.level = CFG$thr$roc$ci, method = "delong"))
    rocrows[[acc]] <- data.frame(disease = dis, cohort = acc, type = "panel",
      n = length(y), auc = as.numeric(r$auc),
      ci_low = ci[1], ci_high = ci[3], genes = length(gp))
    # per-gene AUC too
    for (g in gp) {
      rg <- tryCatch(roc(y, va$X[[g]], quiet = TRUE), error = function(e) NULL)
      if (!is.null(rg)) rocrows[[paste0(acc, g)]] <- data.frame(disease = dis,
        cohort = acc, type = g, n = length(y), auc = as.numeric(rg$auc),
        ci_low = NA, ci_high = NA, genes = 1)
    }
    # ROC plot (panel)
    png(file.path(outdir, sprintf("roc_%s_%s.png", dis, acc)), 950, 950, res = 150)
    plot(r, col = "tomato", lwd = 3, main = sprintf("%s panel ROC — %s (n=%d)", dis, acc, length(y)),
         print.auc = TRUE, print.auc.x = .55, print.auc.y = .15,
         identity.col = "grey70")
    legend("bottomright", bty = "n", cex = .8,
           legend = sprintf("AUC %.3f (95%% CI %.3f-%.3f)", r$auc, ci[1], ci[3]))
    dev.off()
    # DCA
    dca <- dca_df(y, lp, "hub panel")
    ggplot(dca, aes(thr, nb, color = strat)) + geom_line(linewidth = .9) +
      coord_cartesian(ylim = c(-0.05, max(dca$nb) + .02)) +
      scale_color_manual(values = c(`hub panel` = "tomato", `treat all` = "grey50", `treat none` = "black")) +
      labs(title = sprintf("Decision curve — %s / %s", dis, acc),
           x = "threshold probability", y = "net benefit", color = NULL) +
      theme_bw(base_size = 12)
    ggsave(file.path(outdir, sprintf("dca_%s_%s.png", dis, acc)), width = 6, height = 4.6, dpi = 120)
    log_msg(dis, "/", acc, " panel AUC=", signif(r$auc, 3),
            " (", signif(ci[1], 3), "-", signif(ci[3], 3), ") n=", length(y),
            " genes=", length(gp), step = "07")
  }
  if (length(rocrows)) {
    out <- do.call(rbind, rocrows)
    write.csv(out, file.path(outdir, sprintf("07_ROC_%s.csv", dis)), row.names = FALSE)
  }
}

for (dis in active_diseases()) run_dx(dis)
log_msg("07 done.", step = "07")

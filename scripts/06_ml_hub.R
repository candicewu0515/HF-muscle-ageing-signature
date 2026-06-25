# =====================================================================
# 06_ml_hub.R — 6 machine-learning selectors on the concordant shared core
#   features = 26 concordant shared genes; for each disease discovery cohort
#   run LASSO / RF / SVM-RFE / Boruta / XGBoost / glmBoost, vote, then take
#   genes robustly selected in BOTH diseases as the shared diagnostic hubs.
#   outputs:
#     results/06_ml/06_votes_<DIS>.csv
#     results/06_ml/06_hub.csv                 (final shared hubs)
#     results/06_ml/ml_vote_heatmap.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({
  library(glmnet); library(randomForest); library(e1071); library(Boruta)
  library(xgboost); library(caret); library(ggplot2)
  have_mboost <- requireNamespace("mboost", quietly = TRUE)
})
set.seed(CFG$SEED)
qc   <- file.path(CFG$dir$results, "01_qc")
shdir<- file.path(CFG$dir$results, "04_shared")
outdir <- file.path(CFG$dir$results, "06_ml")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

core <- read.csv(file.path(shdir, "04_shared_core.csv"))$gene
log_msg("ML features = ", length(core), " concordant core genes", step = "06")

build_xy <- function(dis) {
  acc <- CFG$datasets[[dis]]$discovery$acc
  d <- readRDS(file.path(qc, sprintf("01_expr_%s_%s.rds", dis, acc)))
  g <- intersect(core, rownames(d$expr))
  X <- t(d$expr[g, , drop = FALSE])
  X <- scale(X)
  y <- factor(d$group, levels = c("control", "case"))
  list(X = X, y = y, genes = g)
}

sel_lasso <- function(X, y) {
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = CFG$thr$ml$nfolds)
  co <- as.matrix(coef(cv, s = "lambda.min"))[-1, 1]
  names(co)[co != 0]
}
sel_rf <- function(X, y) {
  rf <- randomForest(X, y, ntree = CFG$thr$ml$rf_ntree, importance = TRUE)
  imp <- importance(rf, type = 2)[, 1]
  names(imp)[imp > median(imp)]
}
sel_svmrfe <- function(X, y) {
  ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 5, verbose = FALSE)
  fit <- tryCatch(rfe(X, y, sizes = 2:ncol(X), rfeControl = ctrl,
                      method = "svmLinear", metric = "Accuracy"),
                  error = function(e) NULL)
  if (is.null(fit)) {  # fallback: linear SVM weights
    m <- svm(X, y, kernel = "linear"); w <- t(m$coefs) %*% m$SV
    nm <- colnames(X)[order(abs(w[1, ]), decreasing = TRUE)]
    return(head(nm, ceiling(ncol(X) / 2)))
  }
  predictors(fit)
}
sel_boruta <- function(X, y) {
  b <- Boruta(X, y, doTrace = 0, maxRuns = 200)
  f <- TentativeRoughFix(b)
  names(f$finalDecision)[f$finalDecision == "Confirmed"]
}
sel_xgb <- function(X, y) {
  dtrain <- xgb.DMatrix(X, label = as.integer(y) - 1)
  m <- xgb.train(list(objective = "binary:logistic", max_depth = 3, eta = .2),
                 dtrain, nrounds = 60, verbose = 0)
  imp <- xgb.importance(model = m)
  if (is.null(imp) || !nrow(imp)) return(character())
  imp$Feature[imp$Gain > 0]
}
sel_glmboost <- function(X, y) {
  if (!have_mboost) return(NA_character_)   # marks "method unavailable"
  df <- data.frame(y = y, X)
  m <- mboost::glmboost(y ~ ., data = df, family = mboost::Binomial())
  v <- unique(names(coef(m)))
  v <- sub("^X", "", setdiff(v, "(Intercept)"))
  intersect(v, colnames(X))
}

METHODS <- list(LASSO = sel_lasso, RF = sel_rf, `SVM-RFE` = sel_svmrfe,
                Boruta = sel_boruta, XGBoost = sel_xgb, glmBoost = sel_glmboost)

run_ml <- function(dis) {
  dat <- build_xy(dis); X <- dat$X; y <- dat$y
  votes <- matrix(0, length(dat$genes), length(METHODS),
                  dimnames = list(dat$genes, names(METHODS)))
  for (mn in names(METHODS)) {
    s <- tryCatch(METHODS[[mn]](X, y), error = function(e) {
      log_msg("  ", dis, "/", mn, " failed: ", conditionMessage(e), step = "06"); character() })
    if (length(s) == 1 && is.na(s)) { votes[, mn] <- NA; next }   # unavailable
    votes[intersect(s, rownames(votes)), mn] <- 1
  }
  vote_sum <- rowSums(votes, na.rm = TRUE)
  df <- data.frame(gene = rownames(votes), votes, vote = vote_sum,
                   check.names = FALSE)
  df <- df[order(-df$vote), ]
  write.csv(df, file.path(outdir, sprintf("06_votes_%s.csv", dis)), row.names = FALSE)
  n_methods <- sum(!is.na(votes[1, ]))
  log_msg(dis, " ML done (", n_methods, " methods active); >=", CFG$thr$ml$ml_vote,
          " votes: ", sum(vote_sum >= CFG$thr$ml$ml_vote), " genes", step = "06")
  setNames(vote_sum, rownames(votes))
}

vh <- run_ml("HF"); vs <- run_ml("SAR")
g <- intersect(names(vh), names(vs))
tab <- data.frame(gene = g, HF_vote = vh[g], SAR_vote = vs[g])
tab$hub <- tab$HF_vote >= CFG$thr$ml$ml_vote & tab$SAR_vote >= CFG$thr$ml$ml_vote
tab <- tab[order(-tab$hub, -(tab$HF_vote + tab$SAR_vote)), ]
write.csv(tab, file.path(outdir, "06_hub.csv"), row.names = FALSE)
hub <- tab$gene[tab$hub]
log_msg("SHARED HUB genes (>=", CFG$thr$ml$ml_vote, " votes in BOTH): ",
        length(hub), " -> ", paste(hub, collapse = ", "), step = "06")

# vote heatmap
hm <- reshape(tab[, c("gene", "HF_vote", "SAR_vote")], direction = "long",
              varying = c("HF_vote", "SAR_vote"), v.names = "vote",
              times = c("HF", "SAR"), timevar = "disease", idvar = "gene")
hm$gene <- factor(hm$gene, levels = rev(tab$gene))
ggplot(hm, aes(disease, gene, fill = vote)) +
  geom_tile(color = "white") +
  geom_text(aes(label = vote), size = 3) +
  scale_fill_gradient(low = "grey95", high = "tomato", limits = c(0, 6)) +
  labs(title = "ML selection votes (6 methods)", x = NULL, y = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "ml_vote_heatmap.png"), width = 4.6,
       height = 1 + .25 * nrow(tab), dpi = 120, limitsize = FALSE)
log_msg("06 done.", step = "06")

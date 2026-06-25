# =====================================================================
# 12_robustness.R — supplementary evidence addressing peer-review attacks
#   (a) overlap hypergeometric + concordance binomial
#   (b) ML hub label-permutation null (are 11 hubs > chance?)
#   (c) DEG-threshold sensitivity (hub stability)
#   (d) aging-confound: overlap of 26 core with senescence/aging gene sets
#   (e) CLIP length bias: per-hub edges vs 3'UTR length
#   output: results/13_robust/12_*.csv
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(glmnet); library(randomForest); library(xgboost); library(data.table) })
set.seed(CFG$SEED)
R <- CFG$dir$results
out <- file.path(R, "13_robust"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
qc <- file.path(R, "01_qc")

deg <- function(d, lfc = CFG$thr$deg$logFC, p = CFG$thr$deg$padj) {
  t <- read.csv(file.path(R, "02_deg", sprintf("02_DEG_%s.csv", d)))
  s <- subset(t, adj.P.Val < p & abs(logFC) > lfc)
  list(all = t$gene, sig = setNames(sign(s$logFC), s$gene))
}
HF <- deg("HF"); SAR <- deg("SAR")
bg <- intersect(HF$all, SAR$all)
hub <- read.csv(file.path(R, "06_ml", "06_hub.csv")); hub <- hub$gene[hub$hub]

## ---- (a) overlap hypergeometric + concordance binomial ----
ix <- intersect(names(HF$sig), names(SAR$sig)); ov <- length(ix)
con <- sum(HF$sig[ix] == SAR$sig[ix])
p_over <- phyper(ov - 1, length(HF$sig), length(bg) - length(HF$sig), length(SAR$sig), lower.tail = FALSE)
# null concordance prob given marginal up/down rates
pc <- mean(HF$sig > 0) * mean(SAR$sig > 0) + mean(HF$sig < 0) * mean(SAR$sig < 0)
p_con <- binom.test(con, ov, pc, alternative = "greater")$p.value
a <- data.frame(metric = c("overlap", "concordant"),
  observed = c(ov, con), expected = c(length(HF$sig) * length(SAR$sig) / length(bg), ov * pc),
  p = c(p_over, p_con),
  note = c("hypergeometric", sprintf("binomial vs null rate %.2f", pc)))
write.csv(a, file.path(out, "12_overlap_significance.csv"), row.names = FALSE)
log_msg("(a) overlap=", ov, " hyperP=", signif(p_over, 3),
        " | concordant=", con, "/", ov, " binomP=", signif(p_con, 3),
        " (null rate ", round(pc, 2), ")", step = "12")

## ---- (b) ML hub label-permutation null (LASSO+RF+XGB vote>=2, fast) ----
load_xy <- function(dis) {
  d <- readRDS(file.path(qc, sprintf("01_expr_%s_%s.rds", dis, CFG$datasets[[dis]]$discovery$acc)))
  core <- read.csv(file.path(R, "04_shared", "04_shared_core.csv"))$gene
  g <- intersect(core, rownames(d$expr))
  list(X = scale(t(d$expr[g, ])), y = factor(d$group, levels = c("control", "case")), g = g)
}
vote3 <- function(X, y) {
  s <- character()
  cv <- tryCatch(cv.glmnet(X, y, family = "binomial", alpha = 1), error = function(e) NULL)
  if (!is.null(cv)) { co <- as.matrix(coef(cv, s = "lambda.min"))[-1, 1]; s <- c(s, names(co)[co != 0]) }
  rf <- randomForest(X, y, ntree = 300); im <- importance(rf)[, 1]; s <- c(s, names(im)[im > median(im)])
  dt <- xgb.DMatrix(X, label = as.integer(y) - 1)
  m <- xgb.train(list(objective = "binary:logistic", max_depth = 3, eta = .2), dt, nrounds = 40, verbose = 0)
  im2 <- xgb.importance(model = m); if (!is.null(im2) && nrow(im2)) s <- c(s, im2$Feature)
  tb <- table(s); names(tb)[tb >= 2]   # selected by >=2 of 3 fast methods
}
HFd <- load_xy("HF"); SARd <- load_xy("SAR")
obs_hub <- intersect(vote3(HFd$X, HFd$y), vote3(SARd$X, SARd$y))
NPERM <- 200
perm_counts <- integer(NPERM)
for (i in seq_len(NPERM)) {
  set.seed(1000 + i)
  h <- intersect(vote3(HFd$X, sample(HFd$y)), vote3(SARd$X, sample(SARd$y)))
  perm_counts[i] <- length(h)
}
obs_n <- length(obs_hub)
p_ml <- (sum(perm_counts >= obs_n) + 1) / (NPERM + 1)
write.csv(data.frame(observed_hubs = obs_n, perm_mean = mean(perm_counts),
  perm_max = max(perm_counts), p = p_ml, nperm = NPERM),
  file.path(out, "12_ml_permutation.csv"), row.names = FALSE)
log_msg("(b) ML 3-method vote: observed ", obs_n, " shared-selected; perm null mean=",
        round(mean(perm_counts), 1), " p=", signif(p_ml, 3), step = "12")

## ---- (c) DEG-threshold sensitivity (11 hubs survive?) ----
sens <- rbindlist(lapply(c(0.263, 0.378, 0.585, 1.0), function(l) {
  h <- deg("HF", l); s <- deg("SAR", l)
  shared <- intersect(names(h$sig), names(s$sig))
  data.table(logFC = l, n_DEG_HF = length(h$sig), n_DEG_SAR = length(s$sig),
             shared = length(shared),
             hub_retained = sum(hub %in% shared), hub_total = length(hub))
}))
fwrite(sens, file.path(out, "12_threshold_sensitivity.csv"))
log_msg("(c) threshold sensitivity: hub retained in shared set @logFC ",
        paste(sprintf("%.3f:%d/%d", sens$logFC, sens$hub_retained, sens$hub_total), collapse = " "),
        step = "12")

## ---- (d) aging-confound: 26 core vs senescence/aging gene sets ----
aging <- tryCatch({
  suppressMessages(library(msigdbr))
  m <- as.data.frame(msigdbr(species = "Homo sapiens"))
  gs_col <- if ("gs_name" %in% names(m)) "gs_name" else "gs_name"
  sym_col <- if ("gene_symbol" %in% names(m)) "gene_symbol" else "gene_symbol"
  sets <- unique(m[[gs_col]][grepl("SENESCEN|AGING|FRIDMAN|SEN_MAYO", m[[gs_col]], ignore.case = TRUE)])
  uni <- unique(m[[sym_col]])
  agers <- unique(m[[sym_col]][m[[gs_col]] %in% sets])
  core <- read.csv(file.path(R, "04_shared", "04_shared_core.csv"))$gene
  hitc <- intersect(core, agers); hithub <- intersect(hub, agers)
  pc <- phyper(length(hitc) - 1, length(agers), length(uni) - length(agers), length(core), lower.tail = FALSE)
  list(n_sets = length(sets), core_in_aging = length(hitc), core_total = length(core),
       hub_in_aging = length(hithub), genes = paste(hitc, collapse = ";"), p = pc)
}, error = function(e) { log_msg("(d) msigdbr failed: ", conditionMessage(e), step = "12"); NULL })
if (!is.null(aging)) {
  write.csv(data.frame(aging[c("n_sets","core_in_aging","core_total","hub_in_aging","p","genes")]),
            file.path(out, "12_aging_overlap.csv"), row.names = FALSE)
  log_msg("(d) aging/senescence overlap: ", aging$core_in_aging, "/", aging$core_total,
          " core genes (", aging$genes, ") hyperP=", signif(aging$p, 3),
          " | hubs in aging sets=", aging$hub_in_aging, step = "12")
}

## ---- (e) CLIP length bias: per-hub edges vs 3'UTR length ----
clip <- tryCatch({
  suppressMessages({ library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(GenomicFeatures) })
  x <- fread(file.path(R, "10_clip", "mirna_hub_edges.csv"))
  ph <- x[, .(edges = .N, n_miR = uniqueN(miRNA), max_clip = max(clipExpNum)), by = gene]
  eg <- mapIds(org.Hs.eg.db, ph$gene, "ENTREZID", "SYMBOL")
  txl <- transcriptLengths(TxDb.Hsapiens.UCSC.hg38.knownGene, with.utr3_len = TRUE)
  u3 <- tapply(txl$utr3_len, txl$gene_id, max, na.rm = TRUE)
  ph$utr3_len <- as.numeric(u3[as.character(eg)])
  ph$pct_of_edges <- round(100 * ph$edges / sum(ph$edges), 1)
  ph <- ph[order(-edges)]
  ct <- suppressWarnings(cor.test(ph$edges, ph$utr3_len, method = "spearman"))
  fwrite(ph, file.path(out, "12_clip_perhub.csv"))
  list(rho = ct$estimate, p = ct$p.value, ccnd1 = ph[gene == "CCND1", pct_of_edges])
}, error = function(e) { log_msg("(e) CLIP length failed: ", conditionMessage(e), step = "12"); NULL })
if (!is.null(clip))
  log_msg("(e) CLIP per-hub: CCND1=", clip$ccnd1, "% of edges; Spearman(edges,3'UTR len) rho=",
          round(clip$rho, 2), " p=", signif(clip$p, 3), step = "12")
log_msg("12 robustness done.", step = "12")

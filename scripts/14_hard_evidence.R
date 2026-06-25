# =====================================================================
# 14_hard_evidence.R — hard evidence addressing the two core critiques
#   (a) AGE CONFOUND: is the HF signature just chronological aging?
#       - is HF older than controls? (GSE57338 has continuous age)
#       - does the signature discriminate HF AFTER adjusting for age?
#   (b) NON-AGING sub-signature: drop senescence/aging genes, re-test AUC
#   (c) HF SUBTYPE: signature in ischemic vs dilated CM separately; do the
#       two etiologies even share the same direction?
#   output: results/14_hard/14_*.csv
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(pROC); library(limma) })
set.seed(CFG$SEED)
R <- CFG$dir$results
out <- file.path(R, "14_hard"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
qc <- file.path(R, "01_qc")

core <- read.csv(file.path(R, "04_shared", "04_shared_core.csv"))
SGN  <- setNames(sign(core$SAR_logFC), core$gene)
hub  <- read.csv(file.path(R, "06_ml", "06_hub.csv")); hub <- hub$gene[hub$hub]
aging <- strsplit(read.csv(file.path(R, "13_robust", "12_aging_overlap.csv"))$genes, ";")[[1]]
nonaging_hub <- setdiff(hub, aging)
log_msg("hubs=", length(hub), " aging-hubs=", length(intersect(hub, aging)),
        " non-aging hubs=", length(nonaging_hub), " (", paste(nonaging_hub, collapse = ","), ")", step = "14")

sig_score <- function(ex, genes) {
  g <- intersect(genes, intersect(rownames(ex), names(SGN)))
  Z <- t(scale(t(ex[g, , drop = FALSE]))); Z[is.na(Z)] <- 0
  as.numeric(t(Z) %*% SGN[g])
}
aucv <- function(y, s) { r <- roc(y, s, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r, method = "delong")); c(auc = as.numeric(r$auc), lo = ci[1], hi = ci[3]) }

## ---- load GSE57338 with age + etiology ----
d <- readRDS(file.path(qc, "01_expr_HF_GSE57338.rds"))
ex <- d$expr; grp <- d$group
age <- suppressWarnings(as.numeric(as.character(d$meta[["age:ch1"]])))
etio <- as.character(d$meta[["disease status:ch1"]])
y <- as.integer(grp) - 1                                 # case=1

## ---- (a) age confound ----
ok <- !is.na(age)
age_case <- age[ok & y == 1]; age_ctrl <- age[ok & y == 0]
p_age <- wilcox.test(age_case, age_ctrl)$p.value
score8 <- sig_score(ex, hub)
# AUCs: age alone, signature alone, signature residualized on age
a_age   <- aucv(y[ok], age[ok])
a_sig   <- aucv(y[ok], score8[ok])
res_sig <- residuals(lm(score8[ok] ~ age[ok]))           # signature with age removed
a_resid <- aucv(y[ok], res_sig)
# logistic: does score add beyond age?
m_full <- glm(y[ok] ~ score8[ok] + age[ok], family = binomial)
p_score_adj <- coef(summary(m_full))["score8[ok]", "Pr(>|z|)"]
write.csv(data.frame(
  metric = c("median_age_HF", "median_age_ctrl", "age_diff_p",
             "AUC_age_alone", "AUC_signature", "AUC_signature_age_adjusted",
             "signature_p_adjusted_for_age"),
  value = c(median(age_case), median(age_ctrl), p_age,
            a_age["auc"], a_sig["auc"], a_resid["auc"], p_score_adj)),
  file.path(out, "14_age_confound.csv"), row.names = FALSE)
log_msg("(a) HF age med=", median(age_case), " ctrl=", median(age_ctrl), " (p=", signif(p_age, 2),
        ") | AUC age-alone=", signif(a_age["auc"], 3), " sig=", signif(a_sig["auc"], 3),
        " sig|age-adj=", signif(a_resid["auc"], 3), " sig-p-after-age=", signif(p_score_adj, 2), step = "14")

## ---- (b) non-aging sub-signature (external HF + SAR) ----
rows <- list()
for (cf in list(c("HF","GSE5406"), c("SAR","GSE1428"), c("SAR","GSE25941"))) {
  dd <- readRDS(file.path(qc, sprintf("01_expr_%s_%s.rds", cf[1], cf[2])))
  yy <- as.integer(dd$group) - 1
  full <- aucv(yy, sig_score(dd$expr, hub))
  na   <- aucv(yy, sig_score(dd$expr, nonaging_hub))
  rows[[paste(cf, collapse = "_")]] <- data.frame(disease = cf[1], cohort = cf[2],
    AUC_full = full["auc"], AUC_nonaging = na["auc"],
    nonaging_lo = na["lo"], nonaging_hi = na["hi"])
  log_msg("(b) ", cf[1], "/", cf[2], " AUC full(8)=", signif(full["auc"], 3),
          " non-aging(", length(nonaging_hub), ")=", signif(na["auc"], 3),
          " (", signif(na["lo"],3), "-", signif(na["hi"],3), ")", step = "14")
}
write.csv(do.call(rbind, rows), file.path(out, "14_nonaging_auc.csv"), row.names = FALSE)

## ---- (c) HF subtype: ischemic vs dilated ----
sc <- sig_score(ex, hub)
isch <- etio == "ischemic"; dcm <- etio == "idiopathic dilated CMP"; nf <- etio == "non-failing"
ai <- aucv(c(rep(1, sum(isch)), rep(0, sum(nf))), c(sc[isch], sc[nf]))
ad <- aucv(c(rep(1, sum(dcm)),  rep(0, sum(nf))), c(sc[dcm],  sc[nf]))
# do the two etiologies share DEG direction? per-gene logFC ischemic-vs-NF & DCM-vs-NF
sub_logfc <- function(mask) {
  g2 <- factor(ifelse(mask, "case", ifelse(nf, "control", NA)), levels = c("control","case"))
  keep <- !is.na(g2); des <- model.matrix(~ g2[keep])
  fit <- eBayes(lmFit(ex[, keep], des))
  topTable(fit, coef = 2, number = Inf, sort.by = "none")$logFC
}
lf_i <- sub_logfc(isch); lf_d <- sub_logfc(dcm)
rho_etio <- cor(lf_i, lf_d, method = "spearman")
write.csv(data.frame(
  metric = c("AUC_ischemic_vs_NF", "AUC_dilated_vs_NF", "etiology_logFC_spearman"),
  value = c(ai["auc"], ad["auc"], rho_etio)),
  file.path(out, "14_subtype.csv"), row.names = FALSE)
log_msg("(c) signature AUC ischemic-vs-NF=", signif(ai["auc"], 3),
        " dilated-vs-NF=", signif(ad["auc"], 3),
        " | ischemic-vs-dilated transcriptome rho=", signif(rho_etio, 2), step = "14")
log_msg("14 hard-evidence done.", step = "14")

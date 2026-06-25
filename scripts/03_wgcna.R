# =====================================================================
# 03_wgcna.R — WGCNA per disease (discovery cohort)
#   - top-MAD gene filter -> soft power -> blockwiseModules
#   - module-trait correlation (trait = disease status case/control)
#   - lock module(s) most associated with disease
#   output:
#     results/03_wgcna/03_modules_<DIS>.csv       gene -> module
#     results/03_wgcna/03_diseaseModuleGenes_<DIS>.txt
#     results/03_wgcna/heatmap_moduletrait_<DIS>.png
#     results/03_wgcna/softpower_<DIS>.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE)
disableWGCNAThreads()   # macOS foreach %dopar% backend breaks pickSoftThreshold; single-thread is fast at 5k genes
cor <- WGCNA::cor       # avoid stats::cor masking inside blockwiseModules
indir  <- file.path(CFG$dir$results, "01_qc")
outdir <- file.path(CFG$dir$results, "03_wgcna")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
TOPN <- 5000L

run_wgcna <- function(dis) {
  acc <- CFG$datasets[[dis]]$discovery$acc
  d <- readRDS(file.path(indir, sprintf("01_expr_%s_%s.rds", dis, acc)))
  ex <- d$expr; grp <- d$group
  # top-MAD genes
  mad_ <- apply(ex, 1, mad, na.rm = TRUE)
  keep <- names(sort(mad_, decreasing = TRUE))[seq_len(min(TOPN, sum(mad_ > 0)))]
  datExpr <- t(ex[keep, ])                      # samples x genes
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  log_msg(dis, " WGCNA input ", nrow(datExpr), " samples x ", ncol(datExpr), " genes", step = "03")

  # soft threshold
  powers <- CFG$thr$wgcna$powerCandidates
  sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",
                           RsquaredCut = CFG$thr$wgcna$rsqCut, verbose = 0)
  pwr <- sft$powerEstimate
  if (is.na(pwr)) pwr <- if (ncol(datExpr) > 0) ifelse(nrow(datExpr) < 20, 18, 12) else 6
  log_msg(dis, " soft power = ", pwr, " (R2 target ", CFG$thr$wgcna$rsqCut, ")", step = "03")

  png(file.path(outdir, sprintf("softpower_%s.png", dis)), 900, 450, res = 110)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "power", ylab = "scale-free R^2", type = "n", main = paste(dis, "scale-free fit"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, col = "red"); abline(h = CFG$thr$wgcna$rsqCut, col = "blue", lty = 2)
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "power", ylab = "mean connectivity",
       type = "n", main = "mean connectivity")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
  dev.off()

  # modules
  net <- blockwiseModules(datExpr, power = pwr, networkType = "signed",
                          TOMType = "signed", minModuleSize = CFG$thr$wgcna$minModuleSize,
                          mergeCutHeight = CFG$thr$wgcna$mergeCutHeight,
                          deepSplit = CFG$thr$wgcna$deepSplit,
                          numericLabels = TRUE, saveTOMs = FALSE,
                          maxBlockSize = 6000, verbose = 0)
  moduleColors <- labels2colors(net$colors)
  modtab <- data.frame(gene = colnames(datExpr), module = moduleColors)
  write.csv(modtab, file.path(outdir, sprintf("03_modules_%s.csv", dis)), row.names = FALSE)
  log_msg(dis, " modules: ", length(unique(moduleColors)),
          " (", paste(names(sort(table(moduleColors), decreasing = TRUE))[1:3], collapse = ","), " largest)", step = "03")

  # module-trait correlation
  MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
  trait <- data.frame(disease = as.numeric(grp[gsg$goodSamples]) - 1)  # case=1
  rownames(trait) <- rownames(datExpr)
  mtCor <- cor(MEs, trait, use = "p")
  mtP   <- corPvalueStudent(mtCor, nrow(datExpr))

  png(file.path(outdir, sprintf("heatmap_moduletrait_%s.png", dis)), 520, 900, res = 120)
  par(mar = c(4, 9, 3, 2))
  txt <- paste0(signif(mtCor, 2), "\n(", signif(mtP, 1), ")")
  labeledHeatmap(Matrix = mtCor, xLabels = "disease", yLabels = names(MEs),
                 ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
                 textMatrix = txt, setStdMargins = FALSE, cex.text = 0.6,
                 zlim = c(-1, 1), main = paste(dis, "module-trait"))
  dev.off()

  # lock disease module(s): significant & strongest |cor|
  sig <- which(mtP[, 1] < CFG$thr$wgcna$moduleTraitP)
  ord <- sig[order(abs(mtCor[sig, 1]), decreasing = TRUE)]
  topmods <- sub("^ME", "", names(MEs)[ord])
  topmods <- setdiff(topmods, "grey")
  key <- head(topmods, 2)                         # up to 2 strongest disease modules
  dmGenes <- modtab$gene[modtab$module %in% key]
  writeLines(dmGenes, file.path(outdir, sprintf("03_diseaseModuleGenes_%s.txt", dis)))
  log_msg(dis, " disease module(s): ", paste(key, collapse = "+"),
          " | r=", paste(signif(mtCor[ord[seq_along(key)], 1], 2), collapse = ","),
          " | genes=", length(dmGenes), step = "03")
  invisible(list(modules = modtab, key = key, dmGenes = dmGenes))
}

for (dis in active_diseases()) run_wgcna(dis)
log_msg("03 done.", step = "03")

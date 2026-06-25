# =====================================================================
# 01_fetch_qc.R  —  fetch GEO, normalize, probe->gene, group, QC
#   per disease: discovery + validation cohorts
#   outputs:
#     results/01_qc/01_expr_<DIS>_<ACC>.rds   list(expr=gene x sample, meta=)
#     results/01_qc/qc_<DIS>_<ACC>_box.png / _pca.png
# =====================================================================
src <- file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R")
source(src)
suppressMessages({
  library(GEOquery); library(limma); library(ggplot2)
})
options(timeout = 1800)
Sys.setenv(VROOM_CONNECTION_SIZE = 5e6)
outdir <- file.path(CFG$dir$results, "01_qc")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- grouping spec (FILLED after 01a inspection) --------------------
## returns factor with levels c("control","case") or NA to drop sample
## defined in scripts/01_grouping.R for clarity
source(file.path(CFG$root, "scripts", "01_grouping.R"))

## ---- helpers --------------------------------------------------------
# collapse probes -> gene symbol (max-mean probe per gene), log2 if needed
to_gene_matrix <- function(es) {
  ex <- Biobase::exprs(es)
  ex <- ex[rowSums(!is.na(ex)) > 0, , drop = FALSE]
  # log2 transform if looks linear (heuristic from GEO2R), NA-safe
  qx <- as.numeric(quantile(ex, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
  logc <- isTRUE((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0))
  if (logc) { ex[ex < 0] <- NA; ex <- log2(ex + 1) }
  fd <- Biobase::fData(es)
  sym <- NULL; src <- NA
  # 1) direct symbol columns
  for (cand in c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol",
                 "ILMN_Gene", "GeneSymbol", "gene", "GENE")) {
    if (cand %in% colnames(fd)) { sym <- as.character(fd[[cand]]); src <- cand; break }
  }
  # 2) Affy ST gene_assignment: "acc // SYMBOL // desc // ..."  (symbol = 2nd token)
  if (is.null(sym) && "gene_assignment" %in% colnames(fd)) {
    ga <- as.character(fd[["gene_assignment"]])
    sym <- vapply(strsplit(ga, " // ", fixed = TRUE),
                  function(z) if (length(z) >= 2) trimws(z[2]) else NA_character_, "")
    src <- "gene_assignment"
  }
  # 2.5) RefSeq/GenBank accession column -> symbol via org.Hs.eg.db
  if (is.null(sym)) {
    accc <- intersect(c("GB_ACC", "GB_LIST", "GenBank Accession", "RefSeq",
                        "GB_RANGE"), colnames(fd))
    if (length(accc)) {
      suppressMessages(library(org.Hs.eg.db))
      acc0 <- sub("\\..*$", "", as.character(fd[[accc[1]]]))  # strip version
      m <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = acc0, keytype = "REFSEQ",
                                 column = "SYMBOL", multiVals = "first")
      sym <- unname(m[acc0]); src <- paste0(accc[1], "->REFSEQ")
    }
  }
  # 3) RNA-seq / no annotation: rownames may already be symbols or Ensembl
  if (is.null(sym)) {
    rn <- rownames(ex)
    if (mean(grepl("^[A-Za-z][A-Za-z0-9\\-]+$", rn)) > 0.5) { sym <- rn; src <- "rownames" }
    else return(NULL)
  }
  sym <- sub(" ?///.*$", "", sym)       # first symbol of multi-maps
  sym <- trimws(sym)
  sym[sym %in% c("---", "NA", "")] <- NA
  keep <- !is.na(sym) & sym != "" & !is.na(rowMeans(ex, na.rm = TRUE))
  ex <- ex[keep, , drop = FALSE]; sym <- sym[keep]
  # max-mean collapse
  o <- order(rowMeans(ex, na.rm = TRUE), decreasing = TRUE)
  ex <- ex[o, , drop = FALSE]; sym <- sym[o]
  ex <- ex[!duplicated(sym), , drop = FALSE]
  rownames(ex) <- sym[!duplicated(sym)]
  ex <- normalizeBetweenArrays(ex, method = "quantile")
  ex
}

qc_plots <- function(ex, grp, tag) {
  png(file.path(outdir, paste0("qc_", tag, "_box.png")), 1400, 600, res = 110)
  par(mar = c(7, 4, 3, 1))
  boxplot(ex, outline = FALSE, las = 2, cex.axis = .4,
          col = c("steelblue", "tomato")[as.integer(grp)],
          main = paste0(tag, " — normalized expression"))
  dev.off()
  # PCA
  v <- apply(ex, 1, var); top <- names(sort(v, decreasing = TRUE))[1:min(2000, length(v))]
  pc <- prcomp(t(ex[top, ]), scale. = TRUE)
  pv <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], group = grp)
  ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 2.4, alpha = .8) +
    stat_ellipse(level = .9) +
    labs(title = tag, x = paste0("PC1 (", pv[1], "%)"),
         y = paste0("PC2 (", pv[2], "%)")) +
    scale_color_manual(values = c(control = "steelblue", case = "tomato")) +
    theme_bw(base_size = 13)
  ggsave(file.path(outdir, paste0("qc_", tag, "_pca.png")), width = 6, height = 5, dpi = 120)
}

fetch_one <- function(acc, dis) {
  log_msg("fetch ", dis, "/", acc, step = "01")
  es <- getGEO(acc, destdir = CFG$dir$geo, GSEMatrix = TRUE, getGPL = TRUE)
  es <- if (is.list(es)) es[[1]] else es
  pd <- Biobase::pData(es)
  grp <- group_of(acc, pd)                      # from 01_grouping.R
  ex <- to_gene_matrix(es)
  if (is.null(ex)) { log_msg("  !! no symbol annotation for ", acc, step = "01"); return(invisible(NULL)) }
  # align & drop NA-group samples
  keep <- !is.na(grp)
  ex <- ex[, keep, drop = FALSE]; grp <- droplevels(grp[keep])
  tag <- paste0(dis, "_", acc)
  tabs <- table(grp)
  log_msg("  ", acc, " genes=", nrow(ex), " samples=", ncol(ex),
          " | control=", tabs["control"], " case=", tabs["case"], step = "01")
  qc_plots(ex, grp, tag)
  saveRDS(list(expr = ex, group = grp, acc = acc, disease = dis,
               meta = pd[keep, , drop = FALSE]),
          file.path(outdir, paste0("01_expr_", tag, ".rds")))
  invisible(tabs)
}

## ---- driver ---------------------------------------------------------
safe_fetch <- function(acc, dis) tryCatch(fetch_one(acc, dis),
  error = function(e) log_msg("  !! ", dis, "/", acc, " FAILED: ",
                              conditionMessage(e), step = "01"))
for (dis in active_diseases()) {
  d <- CFG$datasets[[dis]]
  safe_fetch(d$discovery$acc, dis)
  for (v in d$validation) safe_fetch(v$acc, dis)
}
log_msg("01 done.", step = "01")

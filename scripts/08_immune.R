# =====================================================================
# 08_immune.R — ssGSEA immune-cell infiltration x hub correlation
#   GSVA::ssgsea on each discovery cohort with canonical immune signatures
#   (open alternative to registration-gated CIBERSORT); correlate immune
#   cell scores with hub gene expression and case/control.
#   outputs per disease:
#     results/08_immune/08_cibersort_<DIS>.csv     (ssGSEA scores)
#     immune_hub_corr_<DIS>.png   immune_group_<DIS>.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(GSVA); library(ggplot2); library(reshape2) })
qc <- file.path(CFG$dir$results, "01_qc")
outdir <- file.path(CFG$dir$results, "08_immune")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

hub <- read.csv(file.path(CFG$dir$results, "06_ml", "06_hub.csv"))
hub <- hub$gene[hub$hub]

# canonical immune-cell marker signatures (Bindea/Charoentong-style, trimmed)
SIG <- list(
  `CD8 T cells`        = c("CD8A","CD8B","GZMA","GZMK","PRF1","EOMES"),
  `CD4 T cells`        = c("CD4","IL7R","CD40LG","FOXP1"),
  `Tregs`              = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"),
  `Th1 cells`          = c("TBX21","IFNG","STAT1","IL12RB2"),
  `Th2 cells`          = c("GATA3","IL4","IL5","IL13"),
  `Tfh cells`          = c("CXCR5","BCL6","ICOS","PDCD1"),
  `B cells`            = c("CD19","MS4A1","CD79A","CD79B","CD22"),
  `Plasma cells`       = c("MZB1","IGHG1","XBP1","PRDM1","SDC1"),
  `NK cells`           = c("NCAM1","KLRD1","KLRF1","NKG7","GNLY"),
  `Monocytes`          = c("CD14","LYZ","S100A8","S100A9","FCN1"),
  `Macrophages M1`     = c("CD68","NOS2","SOCS3","IL1B","CXCL10"),
  `Macrophages M2`     = c("MRC1","CD163","MSR1","C1QA","C1QB"),
  `Dendritic cells`    = c("ITGAX","CD1C","CLEC9A","FCER1A","LAMP3"),
  `pDC`                = c("LILRA4","IL3RA","CLEC4C","GZMB"),
  `Neutrophils`        = c("FCGR3B","CSF3R","CXCR2","S100A12"),
  `Mast cells`         = c("TPSAB1","CPA3","MS4A2","KIT"),
  `Eosinophils`        = c("CCR3","IL5RA","PRG2","SIGLEC8"),
  `Cytotoxic`          = c("GZMB","GZMH","KLRB1","FGFBP2"),
  `Tem`                = c("CCR7","SELL","S1PR1"),
  `MDSC`               = c("ITGAM","ARG1","S100A8","CD33"),
  `Macrophages`        = c("CD68","CD163","CSF1R","AIF1"),
  `T helper cells`     = c("CD4","ANK3","FBLN7")
)

run_imm <- function(dis) {
  acc <- CFG$datasets[[dis]]$discovery$acc
  d <- readRDS(file.path(qc, sprintf("01_expr_%s_%s.rds", dis, acc)))
  ex <- as.matrix(d$expr); grp <- d$group
  par <- ssgseaParam(ex, SIG, minSize = 2)
  es <- gsva(par, verbose = FALSE)                # cells x samples
  write.csv(data.frame(cell = rownames(es), es, check.names = FALSE),
            file.path(outdir, sprintf("08_cibersort_%s.csv", dis)), row.names = FALSE)
  log_msg(dis, " ssGSEA: ", nrow(es), " immune signatures x ", ncol(es), " samples", step = "08")

  # case vs control per cell type (Wilcoxon)
  pv <- apply(es, 1, function(x) tryCatch(wilcox.test(x ~ grp)$p.value, error = function(e) NA))
  gdf <- melt(data.frame(t(es), group = grp, check.names = FALSE),
              id.vars = "group", variable.name = "cell", value.name = "score")
  gdf$cell <- factor(gdf$cell, levels = names(sort(pv)))
  gdf$sig <- ifelse(pv[as.character(gdf$cell)] < 0.05, "*", "")
  ggplot(gdf, aes(cell, score, fill = group)) +
    geom_boxplot(outlier.size = .4, linewidth = .3) +
    scale_fill_manual(values = c(control = "steelblue", case = "tomato")) +
    labs(title = sprintf("%s (%s): immune infiltration case vs control", dis, acc),
         x = NULL, y = "ssGSEA score") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
  ggsave(file.path(outdir, sprintf("immune_group_%s.png", dis)), width = 8, height = 5, dpi = 120)

  # hub x immune correlation heatmap
  hg <- intersect(hub, rownames(ex))
  cm <- cor(t(ex[hg, , drop = FALSE]), t(es), method = "spearman")
  cdf <- melt(cm, varnames = c("hub", "cell"), value.name = "rho")
  ggplot(cdf, aes(cell, hub, fill = rho)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato", limits = c(-1, 1)) +
    labs(title = sprintf("%s: hub gene x immune-cell correlation (Spearman)", dis), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
  ggsave(file.path(outdir, sprintf("immune_hub_corr_%s.png", dis)), width = 9, height = 4.5, dpi = 120)
  log_msg(dis, " immune: ", sum(pv < 0.05, na.rm = TRUE), " cell types differ (p<0.05)", step = "08")
}

for (dis in active_diseases()) run_imm(dis)
log_msg("08 done.", step = "08")

# =====================================================================
# 02_deg.R — limma differential expression (case vs control) per disease
#   input : results/01_qc/01_expr_<DIS>_<discoveryACC>.rds
#   output: results/02_deg/02_DEG_<DIS>.csv  + volcano_<DIS>.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(limma); library(ggplot2) })
indir  <- file.path(CFG$dir$results, "01_qc")
outdir <- file.path(CFG$dir$results, "02_deg")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

run_deg <- function(dis) {
  acc <- CFG$datasets[[dis]]$discovery$acc
  f <- file.path(indir, sprintf("01_expr_%s_%s.rds", dis, acc))
  if (!file.exists(f)) { log_msg("missing ", f, step = "02"); return(NULL) }
  d <- readRDS(f); ex <- d$expr; grp <- d$group
  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  fit <- lmFit(ex, design)
  fit <- eBayes(contrasts.fit(fit, makeContrasts(case - control, levels = design)))
  tt <- topTable(fit, number = Inf, adjust.method = CFG$thr$deg$method)
  tt$gene <- rownames(tt)
  th <- CFG$thr$deg
  tt$sig <- with(tt, ifelse(adj.P.Val < th$padj & logFC >  th$logFC, "up",
                     ifelse(adj.P.Val < th$padj & logFC < -th$logFC, "down", "ns")))
  write.csv(tt, file.path(outdir, sprintf("02_DEG_%s.csv", dis)), row.names = FALSE)
  n <- table(tt$sig)
  log_msg(dis, "/", acc, " DEG up=", n["up"] %||% 0, " down=", n["down"] %||% 0,
          " (|logFC|>", th$logFC, ", padj<", th$padj, ")", step = "02")

  lab <- head(tt[order(tt$adj.P.Val), ][tt$sig != "ns", ], 12)
  ggplot(tt, aes(logFC, -log10(adj.P.Val), color = sig)) +
    geom_point(alpha = .5, size = 1) +
    geom_vline(xintercept = c(-th$logFC, th$logFC), lty = 2, col = "grey50") +
    geom_hline(yintercept = -log10(th$padj), lty = 2, col = "grey50") +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c(up = "tomato", down = "steelblue", ns = "grey80")) +
    labs(title = sprintf("%s (%s): case vs control", dis, acc),
         x = "log2 fold change", y = "-log10 adj.P") +
    theme_bw(base_size = 13)
  ggsave(file.path(outdir, sprintf("volcano_%s.png", dis)), width = 6.5, height = 5.5, dpi = 120)
  invisible(tt)
}

for (dis in active_diseases()) run_deg(dis)
log_msg("02 done.", step = "02")

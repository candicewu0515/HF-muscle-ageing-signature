# =====================================================================
# 04_shared.R — shared gene set across HF & sarcopenia
#   shared      = DEG_HF ∩ DEG_SAR  (symmetric adjP<0.05 & |logFC|>thr)
#   core        = shared genes with SAME direction in both (concordant)
#   module cols = annotation: is gene in the disease-driving WGCNA module?
#   outputs:
#     results/04_shared/04_shared.csv        (all shared + annotation)
#     results/04_shared/04_shared_core.csv   (concordant core -> ML input)
#     results/04_shared/venn_shared.png
#     results/04_shared/direction_consistency.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(ggplot2); library(VennDiagram) })
deg_dir <- file.path(CFG$dir$results, "02_deg")
wg_dir  <- file.path(CFG$dir$results, "03_wgcna")
outdir  <- file.path(CFG$dir$results, "04_shared")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
invisible(file.remove(list.files(outdir, "VennDiagram.*\\.log$", full.names = TRUE)))
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

th <- CFG$thr$deg
get_deg <- function(dis) {
  t <- read.csv(file.path(deg_dir, sprintf("02_DEG_%s.csv", dis)))
  sig <- subset(t, adj.P.Val < th$padj & abs(logFC) > th$logFC)
  list(genes = sig$gene,
       dir = setNames(ifelse(sig$logFC > 0, "up", "down"), sig$gene),
       lfc = setNames(sig$logFC, sig$gene))
}
get_mod <- function(dis) readLines(file.path(wg_dir, sprintf("03_diseaseModuleGenes_%s.txt", dis)))

HF <- get_deg("HF"); SAR <- get_deg("SAR")
modHF <- get_mod("HF"); modSAR <- get_mod("SAR")

shared <- intersect(HF$genes, SAR$genes)
sh <- data.frame(
  gene      = shared,
  HF_logFC  = HF$lfc[shared],  HF_dir  = HF$dir[shared],
  SAR_logFC = SAR$lfc[shared], SAR_dir = SAR$dir[shared],
  in_module_HF  = shared %in% modHF,
  in_module_SAR = shared %in% modSAR
)
sh$concordant <- sh$HF_dir == sh$SAR_dir
sh$mean_absFC <- (abs(sh$HF_logFC) + abs(sh$SAR_logFC)) / 2
sh <- sh[order(-sh$concordant, -sh$mean_absFC), ]
write.csv(sh, file.path(outdir, "04_shared.csv"), row.names = FALSE)

core <- sh[sh$concordant, ]
write.csv(core, file.path(outdir, "04_shared_core.csv"), row.names = FALSE)

log_msg("DEG_HF=", length(HF$genes), " DEG_SAR=", length(SAR$genes),
        " | shared=", nrow(sh), " concordant-core=", nrow(core),
        " (", sum(core$HF_dir == "up"), " up, ", sum(core$HF_dir == "down"), " down)",
        step = "04")
log_msg("core in BOTH disease modules: ",
        sum(core$in_module_HF & core$in_module_SAR), "/", nrow(core), step = "04")

# ---- Venn: DEG_HF vs DEG_SAR ----
venn.diagram(
  x = list(`HF DEG (heart)` = HF$genes, `Sarcopenia DEG (muscle)` = SAR$genes),
  filename = file.path(outdir, "venn_shared.png"),
  imagetype = "png", height = 1500, width = 1800, resolution = 300,
  fill = c("tomato", "steelblue"), alpha = .5, cex = 1.4, cat.cex = 1.0,
  cat.pos = c(-25, 25), cat.dist = c(.05, .05), margin = .12,
  main = "Shared differential genes: HF x sarcopenia")

# ---- direction-consistency scatter ----
ggplot(sh, aes(HF_logFC, SAR_logFC, color = concordant)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_abline(slope = 1, intercept = 0, lty = 3, color = "grey80") +
  geom_point(size = 2.6, alpha = .85) +
  ggrepel::geom_text_repel(data = subset(sh, concordant), aes(label = gene),
                           size = 2.7, max.overlaps = 40, color = "grey20") +
  scale_color_manual(values = c(`TRUE` = "forestgreen", `FALSE` = "grey65"),
                     name = "same direction", labels = c("discordant", "concordant")) +
  labs(title = "Cross-tissue direction consistency of shared genes",
       subtitle = sprintf("%d shared; %d concordant (upper-right + lower-left quadrants)",
                          nrow(sh), nrow(core)),
       x = "HF log2FC (heart / LV)", y = "Sarcopenia log2FC (skeletal muscle)") +
  theme_bw(base_size = 13)
ggsave(file.path(outdir, "direction_consistency.png"), width = 6.8, height = 5.8, dpi = 120)

log_msg("04 done. core (concordant) genes -> 06 ML input", step = "04")

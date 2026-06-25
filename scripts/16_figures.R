# =====================================================================
# 16_figures.R — main-figure assembly (publication).
#   - generates a FOCUSED CLIP core subnetwork (full network -> Supplementary)
#   - Figure 2 (diagnostic): direction / HF ROC / muscle forest / calibration
#   - Figure 3 (mechanism): core miRNA subnetwork / PPI / heart & muscle dotplots
#   output: results/figures/Figure2_diagnostic.png, Figure3_regulation.png,
#           core_mirna_network.png, figure_manifest.md
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(magick); library(data.table); library(igraph) })
R <- CFG$dir$results
out <- file.path(R, "figures"); dir.create(out, showWarnings = FALSE, recursive = TRUE)

## ---- focused CLIP core subnetwork (readable main-figure version) ----
ed <- fread(file.path(R, "10_clip", "mirna_hub_edges.csv"))
core_mir <- c("hsa-miR-15a-5p","hsa-miR-15b-5p","hsa-miR-16-5p","hsa-miR-195-5p","hsa-miR-497-5p",
              "hsa-miR-424-5p","hsa-miR-34a-5p","hsa-miR-26b-5p","hsa-miR-29a-3p","hsa-miR-29b-3p")
core_hub <- c("CCND1","NAMPT","PPM1K","FBLN1","EIF4EBP1","MTFP1")   # incl. length-null survivors
sub <- ed[miRNA %in% core_mir & gene %in% core_hub]
if (nrow(sub) > 0) {
  g <- graph_from_data_frame(sub[, .(from = miRNA, to = gene, w = clipExpNum, val = xref_validated)], directed = TRUE)
  vt <- ifelse(V(g)$name %in% core_hub, "hub", "miRNA")
  png(file.path(out, "core_mirna_network.png"), 2100, 1150, res = 150)
  set.seed(CFG$SEED)
  plot(g, vertex.color = ifelse(vt == "hub", "tomato", "gold"),
       vertex.size = ifelse(vt == "hub", 30, 18), vertex.frame.color = "grey40",
       vertex.label = sub("^hsa-", "", V(g)$name),               # shorten miRNA labels
       vertex.label.cex = ifelse(vt == "hub", 1.0, 0.8), vertex.label.color = "black",
       vertex.label.dist = ifelse(vt == "hub", 0, 0.8),          # offset miRNA labels to reduce overlap
       vertex.label.degree = -pi/2, vertex.label.family = "sans",
       edge.color = ifelse(E(g)$val, "firebrick", "grey60"),
       edge.width = 0.6 + E(g)$w / 12, edge.arrow.size = .4,
       layout = layout_as_bipartite(g, types = vt == "hub"),
       main = "Core CLIP-supported miRNA–hub interactions")
  dev.off()
  log_msg("core subnetwork: ", vcount(g), " nodes, ", ecount(g), " edges", step = "16")
}

lab <- function(path, letter, w = 1000) {
  im <- image_scale(image_read(path), paste0(w))
  im <- image_border(im, "white", "10x44")
  image_annotate(im, letter, size = 52, weight = 700, gravity = "northwest", location = "+4+0", color = "black")
}
row <- function(...) image_append(image_join(...))
col <- function(...) image_append(image_join(...), stack = TRUE)

## Figure 2 — diagnostic (D = calibration, not nomogram)
f2 <- col(
  row(lab(file.path(R, "04_shared", "direction_consistency.png"), "A"),
      lab(file.path(R, "07_diag", "roc_HF_GSE5406.png"), "B")),
  row(lab(file.path(R, "15_meta", "forest_muscle_auc.png"), "C"),
      lab(file.path(R, "07_diag", "calibration_HF.png"), "D"))
)
image_write(f2, file.path(out, "Figure2_diagnostic.png"), density = 150)

## Figure 3 — mechanism: A core network (large, full width) + B/C cell-type dotplots
## (STRING PPI moved to Supplementary for readability)
f3 <- col(
  lab(file.path(out, "core_mirna_network.png"), "A", 2300),
  row(lab(file.path(R, "09_scrna", "dotplot_heart.png"), "B", 1150),
      lab(file.path(R, "09_scrna", "dotplot_muscle.png"), "C", 1150))
)
image_write(f3, file.path(out, "Figure3_regulation.png"), density = 150)

writeLines(c(
  "# Figure manifest — HF x muscle-ageing manuscript",
  "Figure 1 (study design): Figure1_design.svg/.png (vector; refine in Illustrator).",
  "Figure 2 (Diagnostic): A direction-consistency; B HF external ROC (GSE5406, 0.92);",
  "  C muscle 5-cohort meta-AUC forest (DeLong CIs; pooled 0.76); D HF calibration.",
  "Figure 3 (Mechanism): A FOCUSED core miRNA->hub subnetwork (miR-15/16/195/497 -> CCND1/NAMPT/PPM1K/FBLN1;",
  "  edge width proportional to #AGO-CLIP experiments, red = miRTarBase/TarBase-validated);",
  "  B heart hub dotplot; C muscle hub dotplot. Legend should name the key observations:",
  "  NAMPT -> myeloid (both tissues), FBLN1 -> stroma (cardiac fibroblast / tendon).",
  "Supplementary: STRING PPI (05_enrich/ppi_network.png), full CLIP network",
  "  (10_clip/clip_network.png), HF nomogram, DCA,",
  "  volcano x2, WGCNA module-trait x2, Venn, immune ssGSEA, threshold/age/aging/ML/CLIP-null."
), file.path(out, "figure_manifest.md"))
log_msg("16 figures done.", step = "16")

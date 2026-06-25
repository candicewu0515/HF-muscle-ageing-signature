# =====================================================================
# 05_enrich_ppi.R — functional enrichment + STRING PPI on shared genes
#   enrichment on the 59 shared genes (more power); PPI highlights core.
#   outputs:
#     results/05_enrich/05_enrich.csv  (GO-BP/MF/CC, KEGG, Reactome)
#     results/05_enrich/dotplot_<db>.png
#     results/05_enrich/ppi_edges.csv  + ppi_degree.csv  + ppi_network.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
  library(ReactomePA); library(STRINGdb); library(igraph)
})
indir  <- file.path(CFG$dir$results, "04_shared")
outdir <- file.path(CFG$dir$results, "05_enrich")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sh   <- read.csv(file.path(indir, "04_shared.csv"))
core <- read.csv(file.path(indir, "04_shared_core.csv"))$gene
genes <- sh$gene
eg <- bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)
entrez <- eg$ENTREZID
log_msg("enrichment on ", length(genes), " shared genes (", length(entrez), " mapped)", step = "05")

all_res <- list()
save_enrich <- function(obj, tag, showCat = 12) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) { log_msg("  ", tag, ": none", step = "05"); return() }
  df <- as.data.frame(obj); df$db <- tag
  all_res[[tag]] <<- df
  n <- min(showCat, nrow(df))
  p <- dotplot(obj, showCategory = n) + ggtitle(tag) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(outdir, paste0("dotplot_", tag, ".png")), p, width = 7.5,
         height = 1 + .32 * n, dpi = 120, limitsize = FALSE)
  log_msg("  ", tag, ": ", nrow(df), " terms (q<0.05)", step = "05")
}

for (ont in c("BP", "MF", "CC")) {
  ego <- enrichGO(entrez, org.Hs.eg.db, ont = ont, pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2, minGSSize = 3, readable = TRUE)
  save_enrich(ego, paste0("GO_", ont))
}
kk <- tryCatch(enrichKEGG(entrez, organism = "hsa", pvalueCutoff = 0.05),
               error = function(e) NULL)
if (!is.null(kk)) { kk <- setReadable(kk, org.Hs.eg.db, "ENTREZID"); save_enrich(kk, "KEGG") }
rp <- tryCatch(enrichPathway(entrez, pvalueCutoff = 0.05, readable = TRUE),
               error = function(e) NULL)
save_enrich(rp, "Reactome")

if (length(all_res)) {
  out <- do.call(rbind, lapply(all_res, function(d) d[, intersect(
    c("db", "ID", "Description", "GeneRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"),
    names(d))]))
  write.csv(out, file.path(outdir, "05_enrich.csv"), row.names = FALSE)
}

## ---- STRING PPI ----
sdb <- tryCatch(
  STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400,
               input_directory = CFG$dir$geo),
  error = function(e) { log_msg("STRINGdb init failed: ", conditionMessage(e), step = "05"); NULL })
if (!is.null(sdb)) {
  mp <- sdb$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
  hits <- mp$STRING_id
  ia <- tryCatch(sdb$get_interactions(hits), error = function(e) NULL)
  if (!is.null(ia) && nrow(ia) > 0) {
    id2sym <- setNames(mp$gene, mp$STRING_id)
    ed <- data.frame(from = id2sym[ia$from], to = id2sym[ia$to],
                     combined_score = ia$combined_score)
    ed <- unique(ed[!is.na(ed$from) & !is.na(ed$to) & ed$from != ed$to, ])
    write.csv(ed, file.path(outdir, "ppi_edges.csv"), row.names = FALSE)
    g <- graph_from_data_frame(ed, directed = FALSE)
    deg <- sort(degree(g), decreasing = TRUE)
    write.csv(data.frame(gene = names(deg), degree = as.integer(deg)),
              file.path(outdir, "ppi_degree.csv"), row.names = FALSE)
    log_msg("PPI: ", vcount(g), " nodes, ", ecount(g), " edges. top hubs: ",
            paste(head(names(deg), 6), collapse = ","), step = "05")
    V(g)$core <- names(V(g)) %in% core
    png(file.path(outdir, "ppi_network.png"), 1500, 1300, res = 150)
    set.seed(CFG$SEED)
    plot(g, vertex.size = 4 + 2.2 * sqrt(degree(g)),
         vertex.color = ifelse(V(g)$core, "tomato", "lightsteelblue"),
         vertex.frame.color = NA, vertex.label.cex = .6,
         vertex.label.color = "black", edge.color = "grey80",
         layout = layout_with_fr, main = "STRING PPI of shared genes (red = concordant core)")
    dev.off()
  } else log_msg("PPI: no interactions returned", step = "05")
}
log_msg("05 done.", step = "05")

# =====================================================================
# 10_clip_mirna.R  ★ CLIP-anchored miRNA -> hub regulatory network
#   1) ENCORI/starBase miRNA-target for each hub; keep edges with
#      clipExpNum >= N (>= N supporting AGO-CLIP experiments)  <- the anchor
#   2) cross-validate miRNA->hub with multiMiR (miRTarBase/TarBase)
#   3) upstream TFs via ChEA3 (TF -> hub)
#   4) assemble TF–miRNA–hub network, edges annotated with CLIP exp count
#   outputs:
#     results/10_clip/10_clip_net.csv     (all edges + evidence)
#     results/10_clip/mirna_hub_edges.csv  tf_hub_edges.csv
#     results/10_clip/clip_network.png
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(data.table); library(igraph) })
outdir <- file.path(CFG$dir$results, "10_clip")
cdir   <- file.path(CFG$dir$raw, "encori")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
MINCLIP <- CFG$thr$clip$min_ago_clip_exp

hub <- read.csv(file.path(CFG$dir$results, "06_ml", "06_hub.csv"))
hub <- hub$gene[hub$hub]
log_msg("CLIP-miRNA on ", length(hub), " hub genes (clipExpNum >= ", MINCLIP, ")", step = "10")

## ---- 1) ENCORI miRNA-target (cached per gene) ----
fetch_encori <- function(gene) {
  f <- file.path(cdir, paste0(gene, ".tsv"))
  if (!file.exists(f) || file.size(f) < 200) {
    url <- sprintf("https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA=all&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=%s&cellType=all", gene)
    ok <- FALSE
    for (i in 1:4) {
      try({ download.file(url, f, quiet = TRUE, mode = "wb"); ok <- file.size(f) > 200 }, silent = TRUE)
      if (ok) break; Sys.sleep(2 * i)
    }
  }
  dt <- tryCatch(fread(f, skip = "miRNAid", sep = "\t"), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  dt$geneName <- gene
  dt
}
enc <- rbindlist(lapply(hub, function(g) { Sys.sleep(.3); fetch_encori(g) }), fill = TRUE)
log_msg("ENCORI raw edges: ", nrow(enc), step = "10")
enc_f <- enc[clipExpNum >= MINCLIP, .(miRNA = miRNAname, gene = geneName,
              clipExpNum, panCancer = pancancerNum,
              n_program = (PITA > 0) + (RNA22 > 0) + (miRmap > 0) + (microT > 0) +
                          (miRanda > 0) + (PicTar > 0) + (TargetScan > 0))]
enc_f <- unique(enc_f)
log_msg("CLIP-anchored miRNA->hub edges (clipExpNum>=", MINCLIP, "): ", nrow(enc_f),
        " | miRNAs=", length(unique(enc_f$miRNA)), step = "10")

## ---- 2) multiMiR cross-validation (validated: miRTarBase/TarBase) ----
xref <- tryCatch({
  suppressMessages(library(multiMiR))
  res <- get_multimir(org = "hsa", target = hub, table = "validated", summary = FALSE)
  vt <- res@data
  unique(toupper(paste(vt$mature_mirna_id, vt$target_symbol)))
}, error = function(e) { log_msg("multiMiR failed: ", conditionMessage(e), step = "10"); character() })
enc_f$xref_validated <- toupper(paste(enc_f$miRNA, enc_f$gene)) %in% xref
log_msg("miRNA->hub edges also in miRTarBase/TarBase: ",
        sum(enc_f$xref_validated), "/", nrow(enc_f), step = "10")
fwrite(enc_f, file.path(outdir, "mirna_hub_edges.csv"))

## ---- 3) ChEA3 upstream TFs of hub genes ----
tf_edges <- tryCatch({
  payload <- sprintf('{"query_name":"hub","gene_set":[%s]}',
                     paste(sprintf('"%s"', hub), collapse = ","))
  tf_file <- file.path(cdir, "chea3.json")
  for (i in 1:4) {
    rc <- system2("curl", c("-s", "-m", "60", "-X", "POST",
      "-H", shQuote("Content-Type: application/json"),
      "-d", shQuote(payload),
      shQuote("https://maayanlab.cloud/chea3/api/enrich/")), stdout = tf_file)
    if (file.exists(tf_file) && file.size(tf_file) > 100) break; Sys.sleep(2 * i)
  }
  js <- jsonlite::fromJSON(tf_file)
  top <- head(js[["Integrated--meanRank"]], 15)   # top 15 TFs
  rbindlist(lapply(seq_len(nrow(top)), function(i) {
    tf <- top$TF[i]
    tg <- strsplit(top$Overlapping_Genes[i], ",")[[1]]
    data.table(TF = tf, gene = intersect(trimws(tg), hub), rank = i)
  }), fill = TRUE)
}, error = function(e) { log_msg("ChEA3 failed: ", conditionMessage(e), step = "10"); NULL })
if (!is.null(tf_edges) && nrow(tf_edges)) {
  fwrite(tf_edges, file.path(outdir, "tf_hub_edges.csv"))
  log_msg("ChEA3 TF->hub edges: ", nrow(tf_edges), " | TFs=",
          length(unique(tf_edges$TF)), step = "10")
}

## ---- 4) assemble combined network ----
edges <- rbind(
  data.table(from = enc_f$miRNA, to = enc_f$gene, etype = "miRNA->hub",
             weight = enc_f$clipExpNum, validated = enc_f$xref_validated),
  if (!is.null(tf_edges) && nrow(tf_edges))
    data.table(from = tf_edges$TF, to = tf_edges$gene, etype = "TF->hub",
               weight = NA, validated = NA),
  fill = TRUE)
fwrite(edges, file.path(outdir, "10_clip_net.csv"))

# keep network readable: focus miRNAs hitting >=2 hubs OR high CLIP support
mir_deg <- enc_f[, .(nhub = uniqueN(gene), maxclip = max(clipExpNum)), by = miRNA]
keep_mir <- mir_deg[nhub >= 2 | maxclip >= 5, miRNA]
e2 <- edges[(etype == "TF->hub") | (etype == "miRNA->hub" & from %in% keep_mir)]
g <- graph_from_data_frame(e2, directed = TRUE)
vtype <- ifelse(V(g)$name %in% hub, "hub",
         ifelse(V(g)$name %in% enc_f$miRNA, "miRNA", "TF"))
V(g)$type <- vtype
cols <- c(hub = "tomato", miRNA = "gold", TF = "steelblue")
png(file.path(outdir, "clip_network.png"), 1700, 1500, res = 150)
set.seed(CFG$SEED)
plot(g, vertex.color = cols[vtype],
     vertex.size = ifelse(vtype == "hub", 11, 6),
     vertex.frame.color = NA, vertex.label.cex = ifelse(vtype == "hub", .8, .55),
     vertex.label.color = "black",
     edge.color = ifelse(E(g)$etype == "TF->hub", "steelblue", "grey65"),
     edge.arrow.size = .25, edge.width = ifelse(is.na(E(g)$weight), .6, .4 + E(g)$weight / 6),
     layout = layout_with_fr,
     main = sprintf("CLIP-anchored TF–miRNA–hub network (miRNA edges: >=%d AGO-CLIP exp)", MINCLIP))
legend("topright", legend = names(cols), pt.bg = cols, pch = 21, pt.cex = 1.6, bty = "n")
dev.off()

log_msg("network: ", vcount(g), " nodes (", sum(vtype=="hub"), " hub, ",
        sum(vtype=="miRNA"), " miRNA, ", sum(vtype=="TF"), " TF), ",
        ecount(g), " edges", step = "10")
log_msg("10 done.", step = "10")

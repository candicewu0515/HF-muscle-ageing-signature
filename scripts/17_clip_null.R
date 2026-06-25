# =====================================================================
# 17_clip_null.R — 3'UTR-length-matched null for the CLIP-anchored network
#   Tests the reviewer concern: do the hubs attract MORE high-clipExp miRNAs
#   than random genes of similar 3'UTR length, or is it just length/popularity?
#   For each hub, draw length-matched control genes, query ENCORI (clipExpNum>=3),
#   compare miRNA counts. output: results/17_clipnull/17_clip_null.csv + plot
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(data.table); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db); library(GenomicFeatures); library(ggplot2) })
set.seed(CFG$SEED)
R <- CFG$dir$results
out <- file.path(R, "17_clipnull"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
cdir <- file.path(CFG$dir$raw, "encori"); MINCLIP <- CFG$thr$clip$min_ago_clip_exp
hub <- read.csv(file.path(R, "06_ml", "06_hub.csv")); hub <- hub$gene[hub$hub]

## 3'UTR length for all genes
txl <- transcriptLengths(TxDb.Hsapiens.UCSC.hg38.knownGene, with.utr3_len = TRUE)
u3 <- tapply(txl$utr3_len, txl$gene_id, max, na.rm = TRUE)
sym <- mapIds(org.Hs.eg.db, names(u3), "SYMBOL", "ENTREZID")
gene_u3 <- data.table(gene = unname(sym), utr3 = as.numeric(u3))[!is.na(gene) & utr3 > 0]
gene_u3 <- gene_u3[, .(utr3 = max(utr3)), by = gene]

## ENCORI fetch (cached), return # miRNAs with clipExpNum>=MINCLIP
fetch_n <- function(gene) {
  f <- file.path(cdir, paste0(gene, ".tsv"))
  if (!file.exists(f) || file.size(f) < 200) {
    url <- sprintf("https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA=all&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=%s&cellType=all", gene)
    ok <- FALSE
    for (i in 1:3) { try({ download.file(url, f, quiet = TRUE, mode = "wb"); ok <- file.size(f) > 200 }, silent = TRUE); if (ok) break; Sys.sleep(2) }
  }
  dt <- tryCatch(fread(f, skip = "miRNAid", sep = "\t"), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NA_integer_)
  length(unique(dt$miRNAname[dt$clipExpNum >= MINCLIP]))
}

## per hub: observed miRNA count + K length-matched controls
K <- 4
rows <- list()
allgenes <- gene_u3$gene
for (g in hub) {
  hl <- gene_u3[gene == g, utr3]; if (!length(hl)) next
  obs <- fetch_n(g); Sys.sleep(.2)
  pool <- gene_u3[abs(log2(utr3 / hl)) < 0.3 & gene != g & !(gene %in% hub), gene]  # within ~20% length
  ctrl <- sample(pool, min(K, length(pool)))
  cn <- sapply(ctrl, function(c) { n <- fetch_n(c); Sys.sleep(.2); n })
  rows[[g]] <- data.table(hub = g, utr3 = hl, obs_miRNA = obs,
                          ctrl_mean = mean(cn, na.rm = TRUE), ctrl_genes = paste(ctrl, collapse = ";"),
                          ctrl_vals = paste(cn, collapse = ";"))
  log_msg(g, " 3'UTR=", round(hl), " obs=", obs, " ctrl(len-matched) mean=",
          round(mean(cn, na.rm = TRUE), 1), step = "17")
}
res <- rbindlist(rows)
fwrite(res, file.path(out, "17_clip_null.csv"))

## paired test: hubs vs their length-matched controls
ok <- res[!is.na(obs_miRNA) & !is.na(ctrl_mean)]
wt <- tryCatch(wilcox.test(ok$obs_miRNA, ok$ctrl_mean, paired = TRUE), error = function(e) NULL)
enr <- ok[, .(median_hub = median(obs_miRNA), median_ctrl = median(ctrl_mean))]
log_msg("CLIP length-matched null: hub median miRNAs=", enr$median_hub,
        " vs length-matched ctrl=", round(enr$median_ctrl, 1),
        if (!is.null(wt)) paste0(" paired-Wilcoxon p=", signif(wt$p.value, 3)) else "", step = "17")

dl <- melt(ok[, .(hub, observed = obs_miRNA, `length-matched control` = ctrl_mean)],
           id.vars = "hub")
ggplot(dl, aes(variable, value, group = hub)) +
  geom_line(color = "grey75") + geom_point(aes(color = variable), size = 2.5) +
  scale_color_manual(values = c(observed = "tomato", `length-matched control` = "steelblue"), guide = "none") +
  labs(title = "CLIP miRNA edges: hubs vs 3'UTR-length-matched controls",
       subtitle = if (!is.null(wt)) sprintf("paired Wilcoxon p=%.3f", wt$p.value) else "",
       x = NULL, y = "# miRNAs with clipExpNum>=3") +
  theme_bw(base_size = 12)
ggsave(file.path(out, "clip_null.png"), width = 6, height = 5, dpi = 120)
log_msg("17 CLIP null done.", step = "17")

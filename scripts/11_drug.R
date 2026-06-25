# =====================================================================
# 11_drug.R — drug repurposing via DSigDB (Enrichr API) on shared/hub genes
#   reverse-signature candidate drugs whose target sets enrich the shared
#   genes; ranked table + barplot. (Docking/MD NOT run — no vina/gromacs/rdkit
#   in this env, and MD is not feasible in-session; we emit a docking plan
#   pairing top drugs to druggable hub targets for later execution.)
#   outputs:
#     results/11_drug/11_dsigdb.csv     results/11_drug/drug_barplot.png
#     results/11_drug/docking_plan.csv
# =====================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages({ library(jsonlite); library(ggplot2) })
outdir <- file.path(CFG$dir$results, "11_drug")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sh  <- read.csv(file.path(CFG$dir$results, "04_shared", "04_shared.csv"))
hub <- read.csv(file.path(CFG$dir$results, "06_ml", "06_hub.csv"))
hub <- hub$gene[hub$hub]
genes <- unique(sh$gene)                 # use full shared set for enrichment power

## ---- Enrichr: add list, then enrich against DSigDB ----
add_url <- "https://maayanlab.cloud/Enrichr/addList"
enr_url <- "https://maayanlab.cloud/Enrichr/enrich"
tmpdesc <- "HF_SAR_shared"

post_list <- function(genes) {
  body <- paste0("--XX\r\nContent-Disposition: form-data; name=\"list\"\r\n\r\n",
                 paste(genes, collapse = "\n"),
                 "\r\n--XX\r\nContent-Disposition: form-data; name=\"description\"\r\n\r\n",
                 tmpdesc, "\r\n--XX--\r\n")
  f <- tempfile(); writeLines(body, f)
  for (i in 1:4) {
    out <- tryCatch(system2("curl", c("-s", "-m", "40", "-X", "POST",
        "-H", shQuote("Content-Type: multipart/form-data; boundary=XX"),
        "--data-binary", shQuote(paste0("@", f)), shQuote(add_url)), stdout = TRUE),
      error = function(e) "")
    j <- tryCatch(fromJSON(paste(out, collapse = "")), error = function(e) NULL)
    if (!is.null(j$userListId)) return(j$userListId)
    Sys.sleep(2 * i)
  }
  stop("Enrichr addList failed")
}

get_enrich <- function(uid, lib) {
  for (i in 1:4) {
    out <- tryCatch(system2("curl", c("-s", "-m", "40", "-G", shQuote(enr_url),
        "--data-urlencode", shQuote(paste0("userListId=", uid)),
        "--data-urlencode", shQuote(paste0("backgroundType=", lib))), stdout = TRUE),
      error = function(e) "")
    j <- tryCatch(fromJSON(paste(out, collapse = "")), error = function(e) NULL)
    if (!is.null(j[[lib]])) return(j[[lib]])
    Sys.sleep(2 * i)
  }
  NULL
}

uid <- post_list(genes)
log_msg("Enrichr list id=", uid, " (", length(genes), " shared genes)", step = "11")
res <- get_enrich(uid, "DSigDB")
if (is.null(res) || !length(res)) { log_msg("DSigDB returned nothing", step = "11"); quit(save = "no") }

# Enrichr row: [rank, term, pval, zscore, combined, genes, adjp, ...]
df <- do.call(rbind, lapply(res, function(r) data.frame(
  term = r[[2]], pval = as.numeric(r[[3]]), combined = as.numeric(r[[5]]),
  adjp = as.numeric(r[[7]]),
  genes = paste(unlist(r[[6]]), collapse = ";"),
  n = length(unlist(r[[6]])), stringsAsFactors = FALSE)))
df <- df[order(df$adjp, -df$combined), ]
write.csv(df, file.path(outdir, "11_dsigdb.csv"), row.names = FALSE)
log_msg("DSigDB drugs: ", nrow(df), " ; sig(adjP<0.05): ", sum(df$adjp < 0.05),
        " ; top: ", paste(head(df$term, 5), collapse = " | "), step = "11")

top <- head(df[df$adjp < 0.05, ], 20); if (!nrow(top)) top <- head(df, 20)
top$lab <- sub(" .*$", "", top$term)
ggplot(top, aes(reorder(lab, -adjp), -log10(adjp))) +
  geom_col(fill = "mediumpurple") + coord_flip() +
  labs(title = "Candidate drugs (DSigDB) reversing the shared signature",
       x = NULL, y = "-log10 adj.P") + theme_bw(base_size = 11)
ggsave(file.path(outdir, "drug_barplot.png"), width = 7.5, height = 6, dpi = 120)

## ---- docking plan: pair top drugs to druggable hub targets ----
# hubs that hit a drug's gene list become candidate docking targets
plan <- do.call(rbind, lapply(seq_len(min(15, nrow(top))), function(i) {
  hits <- intersect(strsplit(top$genes[i], ";")[[1]], hub)
  if (!length(hits)) return(NULL)
  data.frame(drug = top$lab[i], adjp = top$adjp[i], hub_target = hits)
}))
if (!is.null(plan)) {
  write.csv(plan, file.path(outdir, "docking_plan.csv"), row.names = FALSE)
  log_msg("docking plan: ", nrow(plan), " drug-hub pairs across ",
          length(unique(plan$hub_target)), " hub targets (run vina/MD downstream)",
          step = "11")
}
log_msg("11 done. (DSigDB only; docking/MD deferred — tools absent, MD infeasible in-session)", step = "11")

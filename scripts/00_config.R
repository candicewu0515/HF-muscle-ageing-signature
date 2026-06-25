# =====================================================================
# 00_config.R  —  central config + driver for HF x sarcopenia pipeline
# Sourced by every downstream script (01..11).
#
#   MODE = "single"  -> run one disease alone (set DISEASE)
#   MODE = "cross"   -> run both diseases + shared-set intersection
#
# Usage:
#   Rscript scripts/00_config.R                # prints config, checks env
#   source("scripts/00_config.R")              # from another script
#   MODE=cross Rscript scripts/01_fetch_qc.R   # env override
# =====================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a

## ---- run mode -------------------------------------------------------
CFG <- list()
CFG$MODE    <- Sys.getenv("MODE",    "cross")          # single | cross
CFG$DISEASE <- Sys.getenv("DISEASE", "")               # used when MODE=single: HF | SAR
CFG$SEED    <- 42

## ---- paths ----------------------------------------------------------
# robust root detection: walk up from cwd / script dir until we find scripts/
.find_root <- function() {
  cand <- c(getwd())
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) cand <- c(dirname(normalizePath(f, mustWork = FALSE)),
                           dirname(dirname(normalizePath(f, mustWork = FALSE))), cand)
  for (d in cand) {
    p <- d
    for (i in 1:4) {
      if (dir.exists(file.path(p, "scripts"))) return(normalizePath(p))
      p <- dirname(p)
    }
  }
  getwd()
}
CFG$root <- .find_root()
CFG$dir <- list(
  raw     = file.path(CFG$root, "data", "raw"),
  geo     = file.path(CFG$root, "data", "raw", "geo"),
  results = file.path(CFG$root, "results"),
  logs    = file.path(CFG$root, "logs")
)
for (d in CFG$dir) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- dataset registry (verified 2026-06-25) -------------------------
# group col: "case" label substring(s) vs "control" label substring(s)
# parsed from phenotype in 01; see datasets_manifest.md
CFG$datasets <- list(
  HF = list(
    discovery  = list(acc = "GSE57338", gpl = "GPL11532", tech = "array"),
    validation = list(
      list(acc = "GSE5406", gpl = "GPL96", tech = "array")
      # GSE116250 (RNA-seq) dropped: counts only in supplementary file,
      # series-matrix empty. Re-add via supp processing only if a 2nd
      # external cohort is needed. GSE5406 (N=210) suffices for ROC.
    ),
    scrna = "GSE183852",
    tissue = "left ventricle"
  ),
  SAR = list(
    discovery  = list(acc = "GSE8479", gpl = "GPL2700", tech = "array"),
    validation = list(
      list(acc = "GSE1428",  gpl = "GPL96",  tech = "array"),
      list(acc = "GSE25941", gpl = "GPL570", tech = "array")
    ),
    scrna = "GSE167186",
    tissue = "skeletal muscle"
  )
)

## ---- thresholds -----------------------------------------------------
CFG$thr <- list(
  # 1.3-fold + BH<0.05, SYMMETRIC across both diseases. Rationale: |logFC|>1
  # gives SAR only 9 DEG and a 1-gene shared set; |logFC|>0.378 with adjP<0.05
  # in BOTH independent cohorts + direction-concordance filter yields a robust
  # 59-gene shared / 26-gene concordant core (see results/04_shared).
  deg = list(logFC = 0.378, padj = 0.05, method = "BH"),
  wgcna = list(
    minModuleSize   = 30,
    mergeCutHeight  = 0.25,
    powerCandidates = 1:20,
    rsqCut          = 0.85,      # soft-threshold scale-free fit target
    deepSplit       = 2,
    moduleTraitP    = 0.05       # |module-trait| significance
  ),
  ml = list(
    nfolds   = 10,
    repeats  = 5,
    rf_ntree = 500,
    # 06 uses 6 learners; a gene must be selected by >= ml_vote of them
    ml_vote  = 3
  ),
  clip = list(
    min_ago_clip_exp = 3,        # ENCORI: keep miRNA->hub edges with >= N AGO-CLIP experiments
    require_xref     = TRUE      # also require miRTarBase/TarBase support
  ),
  roc = list(ci = 0.95, boot = 2000)
)

## ---- helpers --------------------------------------------------------
log_msg <- function(..., step = "00") {
  msg <- sprintf("[%s] %s | %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 step, paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(CFG$dir$logs, "run.log"), append = TRUE)
}

# which diseases to run given MODE
active_diseases <- function(cfg = CFG) {
  if (cfg$MODE == "single") {
    stopifnot(cfg$DISEASE %in% c("HF", "SAR"))
    cfg$DISEASE
  } else {
    c("HF", "SAR")
  }
}

set.seed(CFG$SEED)

## ---- self-check when run directly -----------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  if (!interactive()) {
    log_msg("config loaded. MODE=", CFG$MODE,
            "  DISEASE=", ifelse(CFG$DISEASE == "", "(both)", CFG$DISEASE))
    log_msg("active diseases: ", paste(active_diseases(), collapse = ", "))
    log_msg("root: ", CFG$root)
  }
}

# =====================================================================
# run_all.R — master driver. Runs the pipeline in order under one MODE.
#   MODE=cross Rscript scripts/run_all.R            # both diseases + shared
#   MODE=single DISEASE=HF Rscript scripts/run_all.R
#   STEPS=01,02,03 Rscript scripts/run_all.R        # subset of steps
# Heavy/optional steps 09 (scRNA) and 11 (drug/docking) are OFF by default;
# enable with RUN_SCRNA=1 / RUN_DRUG=1.
# =====================================================================
here <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]))
source(file.path(here, "00_config.R"))

ALL <- c(
  "01" = "01_fetch_qc.R",
  "02" = "02_deg.R",
  "03" = "03_wgcna.R",
  "04" = "04_shared.R",
  "05" = "05_enrich_ppi.R",
  "06" = "06_ml_hub.R",
  "07" = "07_diagnostic.R",
  "08" = "08_immune.R",
  "10" = "10_clip_mirna.R"
)
if (Sys.getenv("RUN_SCRNA") == "1") ALL["09"] <- "09_scrna.R"
if (Sys.getenv("RUN_DRUG")  == "1") ALL["11"] <- "11_drug_dock.R"
ALL <- ALL[order(names(ALL))]

want <- Sys.getenv("STEPS", "")
if (nzchar(want)) ALL <- ALL[intersect(strsplit(want, ",")[[1]], names(ALL))]

log_msg("=== run_all START | MODE=", CFG$MODE, " steps=",
        paste(names(ALL), collapse = ","), " ===", step = "RUN")
for (st in names(ALL)) {
  scr <- file.path(here, ALL[st])
  if (!file.exists(scr)) { log_msg("step ", st, " script missing: ", ALL[st], step = "RUN"); next }
  log_msg(">>> step ", st, " : ", ALL[st], step = "RUN")
  t0 <- Sys.time()
  rc <- system2("Rscript", scr, env = paste0("MODE=", CFG$MODE,
                if (CFG$DISEASE != "") paste0(" DISEASE=", CFG$DISEASE) else ""))
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  log_msg("<<< step ", st, ifelse(rc == 0, " OK", paste0(" FAILED rc=", rc)),
          " (", dt, "s)", step = "RUN")
  if (rc != 0) log_msg("   (continuing despite failure; check logs)", step = "RUN")
}
log_msg("=== run_all DONE ===", step = "RUN")

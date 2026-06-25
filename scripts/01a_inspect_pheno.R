# 01a: inspect phenotype of each GEO series to design case/control grouping.
# Prints characteristic columns + unique values so 01 grouping can be exact.
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), "00_config.R"))
suppressMessages(library(GEOquery))
options(timeout = 1200)
Sys.setenv(VROOM_CONNECTION_SIZE = 5e6)

accs <- commandArgs(TRUE)
if (!length(accs)) accs <- c("GSE57338", "GSE8479")

for (acc in accs) {
  cat("\n", strrep("=", 70), "\n", acc, "\n", strrep("=", 70), "\n", sep = "")
  es <- tryCatch(
    getGEO(acc, destdir = CFG$dir$geo, GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) { cat("  DOWNLOAD ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(es)) next
  es <- if (is.list(es)) es[[1]] else es
  pd <- Biobase::pData(es)
  cat("samples:", nrow(pd), " | platform:", annotation(es), "\n")
  # show characteristics + likely grouping columns
  cols <- grep("characteristics|:ch1|title|source", names(pd),
               ignore.case = TRUE, value = TRUE)
  cols <- unique(c("title", cols))
  cols <- cols[cols %in% names(pd)]
  for (cl in cols) {
    u <- unique(as.character(pd[[cl]]))
    if (length(u) <= 25 && length(u) > 1) {
      cat("\n  [", cl, "] ", length(u), " levels:\n", sep = "")
      tb <- sort(table(as.character(pd[[cl]])), decreasing = TRUE)
      for (i in seq_along(tb)) cat(sprintf("      %3d  %s\n", tb[i], names(tb)[i]))
    }
  }
}
cat("\nDONE\n")

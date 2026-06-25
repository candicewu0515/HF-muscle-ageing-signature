# 01_grouping.R — per-accession case/control assignment (verified from 01a)
# returns factor levels c("control","case"); NA = drop sample.
# case  = disease (HF) / aged muscle (sarcopenia proxy)
# control = non-failing LV / young muscle

.fac <- function(x) factor(x, levels = c("control", "case"))

group_of <- function(acc, pd) {
  g <- switch(acc,

    # --- HF: heart failure yes/no ---
    "GSE57338" = {
      v <- as.character(pd[["heart failure:ch1"]])
      ifelse(v == "yes", "case", ifelse(v == "no", "control", NA))
    },
    "GSE5406" = {
      v <- tolower(as.character(pd[["characteristics_ch1"]]))
      ifelse(grepl("heart failure", v), "case",
        ifelse(grepl("normal", v), "control", NA))
    },
    "GSE116250" = {
      v <- tolower(as.character(pd[["disease:ch1"]]))
      ifelse(grepl("non-failing", v), "control",
        ifelse(grepl("cardiomyopathy", v), "case", NA))
    },

    # --- sarcopenia: aged (case) vs young (control) ---
    "GSE8479" = {
      v <- as.character(pd[["Sample Group:ch1"]])
      ifelse(v == "Y", "control", ifelse(v == "O", "case", NA))  # drop OE (post-exercise)
    },
    "GSE1428" = {
      v <- tolower(as.character(pd[["description"]]))
      ifelse(grepl("young", v), "control",
        ifelse(grepl("old", v), "case", NA))
    },
    "GSE25941" = {
      v <- tolower(as.character(pd[["age:ch1"]]))
      ifelse(grepl("young", v), "control",
        ifelse(grepl("old", v), "case", NA))
    },

    stop("no grouping rule for ", acc)
  )
  .fac(g)
}

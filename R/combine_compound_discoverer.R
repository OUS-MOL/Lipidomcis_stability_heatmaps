# =============================================================================
# Combine raw Compound Discoverer "Compounds" exports into the samples-in-rows
# matrix expected by lipid_stability_analysis.R
# =============================================================================
# Authors: Ales Kvasnicka, Barbora Pisklakova
#
# Each input file is a raw Compound Discoverer export: ONE COMPOUND PER ROW
# (Name, m/z, RT, ...) plus one "Norm. Area: <raw filename> (Fxx)" column per
# LC-MS injection. The analysis script needs the opposite layout: one SAMPLE
# per row, one lipid per column. This script:
#
#   - keeps only injections whose file name contains the
#     TEMP_TIMEPOINT_REPLICATE design token (e.g. "-80_0h_1", "RT_48h_3",
#     "-80_3xFT_2") and drops everything else (blanks, pooled QCs, ...)
#   - rewrites each kept injection to a clean id "<SAMPLE_PREFIX>_<temp>_<time>_<rep>"
#   - tags every lipid name with its file tag, e.g. "[POS]" / "[NEG]", so
#     annotations that overlap between ionisation modes are never silently merged
#   - de-duplicates repeated Name values within one file (isomers / several
#     features sharing one putative ID) with " (1)", " (2)", ... suffixes,
#     added AFTER the class text so that the lipid-class parser in the main
#     script (text before the first "(") is unaffected
#   - writes the samples-in-rows matrix to OUTPUT_FILE
#
# Run from the repository root:
#   Rscript R/combine_compound_discoverer.R
# =============================================================================

# ---- Configuration (edit here) ------------------------------------------------
INPUT_FILES <- c(                       # tag = path; one entry per ionisation mode
  POS = "data/example/example_POS.xlsx",
  NEG = "data/example/example_NEG.xlsx"
)
OUTPUT_FILE   <- "data/example/example_combined_POS_NEG.xlsx"
SAMPLE_PREFIX <- "EXAMPLE"              # prefix of the clean sample ids

# TEMP_TIMEPOINT_REPLICATE pattern; adjust if your design uses other labels
DESIGN_PAT <- "(-80|-20|\\+4|RT)_(0h|2h|16h|24h|48h|140h|2xFT|3xFT|4xFT|5xFT)_(\\d+)"
# -------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr); library(openxlsx)
})

load_polarity <- function(file, tag) {
  d <- read_excel(file, sheet = 1)

  area_cols <- grep("^Norm\\. Area: ", names(d), value = TRUE)
  raw_name  <- str_remove(area_cols, "^Norm\\. Area: ")
  m         <- str_match(raw_name, DESIGN_PAT)
  keep      <- !is.na(m[, 1])

  message(sprintf("[%s] %d/%d injections match the design; dropping: %s",
                  tag, sum(keep), length(area_cols),
                  paste(raw_name[!keep], collapse = ", ")))

  area_cols <- area_cols[keep]
  clean_id  <- paste0(SAMPLE_PREFIX, "_", m[keep, 2], "_", m[keep, 3], "_", m[keep, 4])
  stopifnot(!anyDuplicated(clean_id))

  sub <- d %>% select(Name, all_of(area_cols))
  names(sub)[-1] <- clean_id

  dup_idx <- which(duplicated(sub$Name) | duplicated(sub$Name, fromLast = TRUE))
  if (length(dup_idx) > 0) {
    occurrence <- ave(seq_along(sub$Name), sub$Name, FUN = seq_along)
    sub$Name[dup_idx] <- paste0(sub$Name[dup_idx], " (", occurrence[dup_idx], ")")
    message(sprintf("[%s] de-duplicated %d repeated Name value(s)", tag, length(dup_idx)))
  }

  sub$Name <- paste0(sub$Name, " [", tag, "]")

  sub %>%
    pivot_longer(-Name, names_to = "sample_id", values_to = "area") %>%
    pivot_wider(names_from = Name, values_from = area)
}

wide_list <- Map(load_polarity, INPUT_FILES, names(INPUT_FILES))

ids <- lapply(wide_list, function(w) w$sample_id)
stopifnot(all(vapply(ids, function(x) setequal(x, ids[[1]]), logical(1))))

combined <- Reduce(function(a, b) left_join(a, b, by = "sample_id"), wide_list)

stopifnot(!any(is.na(combined)))
stopifnot(!any(combined[, -1] == 0))

message(sprintf("\nCombined matrix: %d samples x %d lipid features (%s)",
                nrow(combined), ncol(combined) - 1,
                paste(sprintf("%d %s", vapply(wide_list, ncol, integer(1)) - 1,
                              names(wide_list)), collapse = " + ")))

names(combined)[1] <- "Sample"
write.xlsx(combined, OUTPUT_FILE, overwrite = TRUE)
message("Wrote: ", OUTPUT_FILE)

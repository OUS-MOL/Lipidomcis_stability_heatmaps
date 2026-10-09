# =============================================================================
# Generate FICTIONAL example data for the lipid stability pipeline
# =============================================================================
# Authors: Ales Kvasnicka, Barbora Pisklakova
#
# Writes two synthetic "Compound Discoverer"-style exports (one compound per
# row, one "Norm. Area: <raw file> (Fxx)" column per injection) to
# data/example/:
#
#   example_POS.xlsx   positive ionisation mode
#   example_NEG.xlsx   negative ionisation mode
#
# EVERYTHING HERE IS INVENTED. Lipid names follow common shorthand
# nomenclature but the species list, m/z, retention times and peak areas are
# random numbers with a built-in, artificial stability signal. Nothing in this
# file comes from, or can be traced back to, a real measurement.
#
# Simulated signal (log2 scale, relative to the -80 C / 0 h baseline):
#   * oxidised lipids (name contains CHO / COOH) rise with storage time,
#     faster at RT than at +4 C, barely at -20 C; they also rise with
#     freeze-thaw cycles
#   * lysophospholipids (LPC / LPE) rise with storage time and with
#     freeze-thaw cycles
#   * polyunsaturated PC / PE species fall with storage time and with
#     freeze-thaw cycles
#   * polyunsaturated TG species are stable over time but fall with
#     freeze-thaw cycles
#   * ~8 % of the remaining lipids get a small random time/temperature drift
#   * everything else is flat; replicate noise is ~12 % (log2 sd 0.12)
# The effect sizes are deliberately large so that every figure and statistical
# branch of the pipeline has something to show; they do not describe real
# lipid behaviour.
#
# Besides the 60 design injections each file also contains blank, pooled-QC
# and bare-QC injections, so the filtering done by
# R/combine_compound_discoverer.R has something to remove. A few species are
# listed twice in POS (isomer-like duplicates) and several PC species appear
# in both polarities, so the de-duplication and [POS]/[NEG] tagging have
# something to do as well.
#
# Run from the repository root:
#   Rscript scripts/generate_example_data.R
# =============================================================================

suppressPackageStartupMessages(library(openxlsx))

set.seed(20261009)

SAMPLE_PREFIX <- "EXAMPLE"
OUT_DIR       <- file.path("data", "example")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Experimental design --------------------------------------------------
temps <- c("-20", "+4", "RT")
times <- c("2h", "16h", "24h", "48h", "140h")
fts   <- c("2xFT", "3xFT", "4xFT", "5xFT")

design <- rbind(
  data.frame(temp = "-80", time = "0h"),
  expand.grid(temp = temps, time = times, stringsAsFactors = FALSE),
  data.frame(temp = "-80", time = fts)
)
design  <- design[, c("temp", "time")]
samples <- do.call(rbind, lapply(1:3, function(r) cbind(design, rep = r)))
samples$label  <- sprintf("%s_%s_%s_%d", SAMPLE_PREFIX, samples$temp, samples$time, samples$rep)
samples$hours  <- ifelse(grepl("h$", samples$time),
                         suppressWarnings(as.numeric(sub("h$", "", samples$time))), 0)
samples$cycles <- ifelse(grepl("xFT$", samples$time),
                         suppressWarnings(as.numeric(sub("xFT$", "", samples$time))), 0)
samples$kind   <- "design"

extras <- data.frame(
  label  = c("Blank_processing_1",
             "PQC_1_1", "PQC_1_2", "PQC_1_3", "PQC_2_1", "PQC_2_2",
             paste0(SAMPLE_PREFIX, "_", 1:3)),
  temp   = "-80", time = "0h", rep = NA_integer_,
  hours  = 0, cycles = 0,
  kind   = c("blank", rep("qc", 8)),
  stringsAsFactors = FALSE
)
injections <- rbind(samples[, names(extras)], extras)
stopifnot(sum(injections$kind == "design") == 60)

# ---- 2. Fictional lipid species ----------------------------------------------
fa <- c("14:0", "16:0", "16:1", "18:0", "18:1", "18:2",
        "20:3", "20:4", "22:5", "22:6")

sample_pairs <- function(cls, n) {
  grid <- t(combn(fa, 2)); idx <- sample(nrow(grid), n)
  sprintf("%s(%s_%s)", cls, grid[idx, 1], grid[idx, 2])
}
sample_triples <- function(cls, n) {
  grid <- t(combn(fa, 3)); idx <- sample(nrow(grid), n)
  sprintf("%s(%s_%s_%s)", cls, grid[idx, 1], grid[idx, 2], grid[idx, 3])
}

pc <- sample_pairs("PC", 14)

pos_names <- c(
  sample_triples("TG", 24),
  sprintf("TG(O-%s)", c("50:1", "50:2", "52:2", "52:3", "54:3", "54:4", "56:5", "58:7")),
  sample_pairs("DG", 5),
  sprintf("ChE(%s)", c("16:0", "18:1", "18:2", "20:4")),
  sprintf("ZyE(%s)", c("16:0", "18:1")),
  sprintf("AcCa(%s)", c("12:0", "14:0", "16:0", "18:1", "18:2")),
  c("SM(d18:1_16:0)", "SM(d18:1_18:0)", "SM(d18:1_24:0)",
    "SM(d18:1_24:1)", "SM(t18:0_16:0)", "SM(d18:2_16:0)"),
  c("Cer(d18:1_16:0)", "Cer(d18:1_24:0)", "Cer(d18:1_24:1)",
    "Cer(t18:0_24:0)", "Cer(m18:0_22:0)"),
  sprintf("LPC(%s)", c("14:0", "16:0", "16:1", "18:0", "18:1", "18:2", "20:4")),
  sprintf("LPC(O-%s)", c("16:0", "18:0", "18:1")),
  pc,
  sprintf("PC(O-%s)", c("32:0", "32:1", "34:1", "34:2", "36:2", "36:4", "38:5", "38:6")),
  sprintf("PC(P-%s)", c("34:1", "36:4", "38:4", "38:6")),
  c("PC(16:0_9:0(CHO))", "PC(16:0_9:0(COOH))", "PC(18:0_9:0(CHO))", "PC(18:1_9:0(COOH))"),
  c("TG(16:0_18:1_9:0(CHO))", "TG(16:0_18:1_9:0(COOH))",
    "TG(18:1_18:1_9:0(CHO))", "TG(18:0_18:1_9:0(COOH))")
)
stopifnot(!anyDuplicated(pos_names))
# two isomer-like duplicates: same putative ID, separate feature
pos_names <- c(pos_names, pos_names[1], pc[1])

neg_names <- c(
  sample_pairs("PE", 8),
  sprintf("PE(O-%s)", c("36:4", "38:5", "38:6", "40:6")),
  sprintf("PE(P-%s)", c("36:4", "38:4", "38:6", "40:6")),
  sprintf("LPE(%s)", c("16:0", "18:0", "18:1", "20:4")),
  c("CerP(d18:1_16:0)", "CerP(d18:1_24:0)"),
  sample_pairs("PI", 3),
  c("PE(18:0_9:0(CHO))", "PE(18:0_9:0(COOH))"),
  pc[1:6]                                   # annotated in both polarities
)
stopifnot(!anyDuplicated(neg_names))

# ---- 3. Signal model ---------------------------------------------------------
TEMP_FACTOR <- c("-80" = 0, "-20" = 0.08, "+4" = 0.35, "RT" = 1)
time_kernel <- function(h) log2(1 + h / 12) / log2(1 + 140 / 12)

lipid_params <- function(nm) {
  is_ox   <- grepl("CHO|COOH", nm)
  is_lyso <- grepl("^LP[CE]\\(", nm)
  is_pufa <- grepl("^(PC|PE)\\(.*(20:4|22:5|22:6)", nm) |
             grepl("^(PC|PE)\\((O-|P-)[0-9]+:[4-6]\\)", nm)
  is_tg_pufa <- grepl("^TG\\(.*(20:4|22:5|22:6)", nm)
  if (is_ox)      return(list(A =  runif(1, 1.8, 2.6), ft =  0.30))
  if (is_lyso)    return(list(A =  runif(1, 1.0, 1.8), ft =  0.25))
  if (is_pufa)    return(list(A = -runif(1, 0.9, 1.5), ft = -0.25))
  if (is_tg_pufa) return(list(A = 0,                   ft = -0.25))
  if (runif(1) < 0.08)
    return(list(A = sample(c(-1, 1), 1) * runif(1, 0.6, 1.0), ft = 0))
  list(A = 0, ft = 0)
}

simulate_areas <- function(lipid_names, inj) {
  n   <- length(lipid_names)
  par <- lapply(lipid_names, lipid_params)
  A   <- vapply(par, function(p) p$A,  numeric(1))
  ft  <- vapply(par, function(p) p$ft, numeric(1))
  mu  <- rnorm(n, mean = 6.8, sd = 0.6)                       # log10 baseline area
  tf  <- unname(TEMP_FACTOR[inj$temp])
  eff <- outer(tf * time_kernel(inj$hours), A) + outer(inj$cycles, ft)
  noise <- matrix(rnorm(length(eff), 0, 0.12), nrow(eff)) +   # per-lipid replicate noise
           rnorm(nrow(eff), 0, 0.05)                          # per-injection offset
  areas <- sweep(2^(eff + noise), 2, 10^mu, "*")
  areas[inj$kind == "blank", ] <- areas[inj$kind == "blank", ] * 0.02
  t(areas)                                                    # lipids x injections
}

write_export <- function(lipid_names, tag, file) {
  inj_order <- injections[sample(nrow(injections)), ]
  areas     <- simulate_areas(lipid_names, inj_order)
  fnums     <- sample(1:90, nrow(inj_order))
  colnames(areas) <- sprintf("Norm. Area: %s_fullMS_%s.raw (F%d)",
                             inj_order$label, tag, fnums)
  out <- data.frame(
    Name = lipid_names,
    `m/z` = round(runif(length(lipid_names), 300, 1100), 4),
    `RT [min]` = round(runif(length(lipid_names), 2, 22), 2),
    check.names = FALSE
  )
  out <- cbind(out, as.data.frame(areas, check.names = FALSE))
  write.xlsx(out, file, sheetName = "Compounds", overwrite = TRUE)
  message(sprintf("Wrote %s: %d compounds x %d injections",
                  file, nrow(out), nrow(inj_order)))
}

write_export(pos_names, "POS", file.path(OUT_DIR, "example_POS.xlsx"))
write_export(neg_names, "NEG", file.path(OUT_DIR, "example_NEG.xlsx"))

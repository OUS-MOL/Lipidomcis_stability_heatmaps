# =============================================================================
# Lipid stability analysis: log2 fold change, Kruskal-Wallis / Dunn's tests and
# multi-panel ComplexHeatmap figures
# =============================================================================
# Authors: Ales Kvasnicka, Barbora Pisklakova
#
# Samples stored under different temperature / freeze-thaw conditions are
# compared with a -80 C / 0 h baseline. For every lipid (and for every lipid
# class) the script computes log2 fold changes, tests for differences and draws
# four side-by-side heatmap panels: -20 C | +4 C | RT | freeze-thaw.
#
# Input layout: ONE SAMPLE PER ROW, one lipid per column.
#   column 1     sample name containing TEMP_TIMEPOINT_REPLICATE
#                (e.g. "EXAMPLE_-80_0h_1"; any prefix / suffix is ignored)
#   columns 2..  peak areas (no missing values, no zeros)
# Raw Compound Discoverer exports must be reshaped first, see
# R/combine_compound_discoverer.R.
#
# Usage (from the repository root):
#   Rscript R/lipid_stability_analysis.R                       # example data -> output/
#   Rscript R/lipid_stability_analysis.R my_data.xlsx my_out   # own data
# or open the .Rproj in RStudio and source() this file.
# =============================================================================

# ---- Input / output -------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
INPUT_FILE <- if (length(args) >= 1) args[1] else "data/example/example_combined_POS_NEG.xlsx"
OUTPUT_DIR <- if (length(args) >= 2) args[2] else "output"

# ---- Experimental design --------------------------------------------------------
BASELINE_TEMP <- "-80"
BASELINE_TIME <- "0h"
TEMP_ORDER    <- c("-20", "+4", "RT")
TIME_ORDER    <- c("2h", "16h", "24h", "48h", "140h")
FT_ORDER      <- c("2xFT", "3xFT", "4xFT", "5xFT")
N_REPLICATES  <- 3          # expected replicates per condition (sanity check)

# ---- Comparison groups ----------------------------------------------------------
# Kruskal-Wallis and Dunn's test run separately per group; the baseline is
# added to every group automatically.
COMPARISON_GROUPS <- list(
  Timepoint   = as.vector(outer(TEMP_ORDER, TIME_ORDER, paste, sep = "_")),
  Freeze_Thaw = paste0(BASELINE_TEMP, "_", FT_ORDER)
)

# ---- Heatmap panels (display layout only) -----------------------------------
HEATMAP_PANELS <- list(
  "-20" = paste0("-20_", TIME_ORDER),
  "+4"  = paste0("+4_",  TIME_ORDER),
  "RT"  = paste0("RT_",  TIME_ORDER),
  "FT"  = paste0(BASELINE_TEMP, "_", FT_ORDER)
)
# Text shown above each panel when it differs from the key in HEATMAP_PANELS
DEG <- intToUtf8(176L)   # degree sign, written this way to keep the source ASCII-only
PANEL_DISPLAY_LABELS <- c("-20" = paste0("-20", DEG, "C"), "+4" = paste0("+4", DEG, "C"))

# ---- Lipid classes --------------------------------------------------------------
# Class = text before the first "(", plus " O-" / " P-" for ether lipids and an
# "ox" prefix when the name contains CHO or COOH. CLASS_MERGE_MAP then relabels
# classes after that parsing (use c() for no merging).
CLASS_MERGE_MAP <- c(
  "ZyE"   = "ChE/ZyE",  "ChE"   = "ChE/ZyE",
  "CerP"  = "Cer",
  "PC O-" = "PC O-/P-", "PC P-" = "PC O-/P-",
  "PE O-" = "PE O-/P-", "PE P-" = "PE O-/P-",
  "TG O-" = "TG"
)
EXCLUDE_OXIDIZED <- TRUE    # drop all classes starting with "ox" before any statistics

# ---- Misc -----------------------------------------------------------------------
SAMPLE_NAME_REGEX <- NULL   # NULL = build automatically from the design above
HEATMAP_LIMITS    <- c(-2, 2)   # log2FC colour scale (values are clipped)
SIG_THRESHOLD     <- 0.05       # BH-adjusted p-value cut-off
TOP_N_PER_CLASS   <- 10         # top-N class summary: lipids kept per class ...
TOP_N_RANK_GROUP  <- "Timepoint" # ... ranked by this group's BH-adjusted KW p-value

set.seed(42)

# =============================================================================
# PACKAGES
# =============================================================================
if (.Platform$OS.type == "windows")
  options(install.packages.compile.from.source = "never")

cran_pkgs <- c("readxl", "openxlsx", "dplyr", "tidyr", "tibble",
               "stringr", "purrr", "circlize", "FSA")
for (p in cran_pkgs)
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap", ask = FALSE)
}

suppressPackageStartupMessages({
  library(readxl);   library(openxlsx); library(dplyr);  library(tidyr)
  library(tibble);   library(stringr);  library(purrr);  library(ComplexHeatmap)
  library(circlize); library(grid);     library(FSA)
})

# =============================================================================
# 1. PARSING HELPERS
# =============================================================================

build_sample_regex <- function() {
  if (!is.null(SAMPLE_NAME_REGEX)) return(SAMPLE_NAME_REGEX)

  esc <- function(x) str_replace_all(x, fixed("+"), "\\+")

  # longest tokens first so that "140h" is tried before "14h", "2xFT" before "2h"
  all_temps <- c(BASELINE_TEMP, TEMP_ORDER)
  all_temps <- all_temps[order(nchar(all_temps), decreasing = TRUE)]
  all_times <- c(BASELINE_TIME, TIME_ORDER, FT_ORDER)
  all_times <- all_times[order(nchar(all_times), decreasing = TRUE)]

  temp_pat <- paste(sapply(all_temps, esc), collapse = "|")
  time_pat <- paste(all_times, collapse = "|")

  paste0("(", temp_pat, ")_(", time_pat, ")_(\\d+)")
}

parse_sample_names <- function(sample_names) {
  regex <- build_sample_regex()
  m     <- str_match(sample_names, regex)

  if (any(is.na(m[, 1]))) {
    bad <- sample_names[is.na(m[, 1])]
    stop(sprintf("Could not parse %d sample name(s) with the current regex.\n",
                 length(bad)),
         "First unparsed name: ", bad[1], "\n",
         "Check TEMP_ORDER / TIME_ORDER / FT_ORDER / BASELINE_TEMP match your names,\n",
         "or override SAMPLE_NAME_REGEX with a custom pattern.")
  }

  tibble(
    sample_id       = sample_names,
    temperature     = m[, 2],
    timepoint_or_FT = m[, 3],
    replicate       = as.integer(m[, 4]),
    condition       = paste0(m[, 2], "_", m[, 3])
  )
}

# extract_lipid_class(): base class -> ether suffix -> oxidised prefix -> merge map
#   "PC(16:0_18:1)"             -> "PC"
#   "PC(O-34:1)"                -> "PC O-"   (then merged to "PC O-/P-")
#   "PE(18:0_9:0(CHO))"         -> "oxPE"
#   "Cer(d18:1_24:0)" / "Cer(t18:0_24:0)" -> both "Cer" (d/t/m sit inside the brackets)
extract_lipid_class <- function(lipid_names) {
  vapply(lipid_names, function(nm) {

    base <- str_trim(str_extract(nm, "^[^(]+"))
    if (is.na(base) || base == "") return(NA_character_)

    ether_suffix <- ""
    m_pl <- str_match(nm, "^[^(]+\\(([OP])-")
    if (!is.na(m_pl[1, 1])) ether_suffix <- paste0(" ", m_pl[1, 2], "-")

    is_oxidized <- str_detect(nm, "CHO|COOH")

    cls <- paste0(base, ether_suffix)
    if (!is.na(is_oxidized) && is_oxidized) cls <- paste0("ox", cls)

    if (cls %in% names(CLASS_MERGE_MAP)) cls <- CLASS_MERGE_MAP[[cls]]

    cls
  }, character(1L), USE.NAMES = FALSE)
}

# Append a/b/c/... to duplicated display labels so every heatmap row is unique.
deduplicate_labels <- function(labels) {
  result  <- labels
  dup_nms <- names(which(table(labels) > 1))
  for (nm in dup_nms) {
    idx         <- which(labels == nm)
    result[idx] <- paste0(nm, letters[seq_along(idx)])
  }
  result
}

# =============================================================================
# 2. DATA LOADING
# =============================================================================

load_data <- function(input_file, exclude_oxidized = FALSE) {
  message("Reading input file: ", input_file)
  raw <- read_excel(input_file)

  sample_col  <- names(raw)[1]
  lipid_names <- names(raw)[-1]

  lipid_class_map <- tibble(
    lipid       = lipid_names,
    lipid_class = extract_lipid_class(lipid_names)
  )

  if (exclude_oxidized) {
    dropped <- lipid_class_map %>% filter(str_starts(lipid_class, "ox")) %>% pull(lipid)
    if (length(dropped) > 0) {
      message(sprintf("Excluding %d oxidized lipid(s) (classes: %s)",
                      length(dropped),
                      paste(sort(unique(lipid_class_map$lipid_class[lipid_class_map$lipid %in% dropped])),
                            collapse = ", ")))
      raw             <- raw %>% select(-all_of(dropped))
      lipid_names     <- setdiff(lipid_names, dropped)
      lipid_class_map <- lipid_class_map %>% filter(!lipid %in% dropped)
    }
  }
  lipid_order <- lipid_names

  meta <- parse_sample_names(raw[[sample_col]])

  lipid_display_labels <- deduplicate_labels(lipid_names)
  n_dups <- sum(lipid_display_labels != lipid_names)
  if (n_dups > 0)
    message(sprintf("Deduplicated %d feature name(s) with a/b/c suffixes", n_dups))

  long_df <- raw %>%
    rename(sample_id = 1) %>%
    pivot_longer(cols = -sample_id, names_to = "lipid", values_to = "raw_value") %>%
    left_join(meta, by = "sample_id") %>%
    left_join(lipid_class_map, by = "lipid") %>%
    mutate(
      lipid       = factor(lipid, levels = lipid_order),
      lipid_class = factor(lipid_class, levels = unique(lipid_class_map$lipid_class)),
      log10_value = log10(raw_value)
    )

  if (any(is.na(long_df$raw_value))) stop("NA values detected in peak area data.")
  if (any(long_df$raw_value == 0))   stop("Zero values detected in peak area data.")

  message(sprintf("Loaded %d samples, %d lipids, %d lipid classes, %d conditions",
                  n_distinct(long_df$sample_id), n_distinct(long_df$lipid),
                  n_distinct(long_df$lipid_class, na.rm = TRUE),
                  n_distinct(long_df$condition)))
  cls_tbl <- lipid_class_map %>% count(lipid_class, sort = TRUE)
  message("Class breakdown:")
  for (i in seq_len(nrow(cls_tbl)))
    message(sprintf("  %-20s  n = %d", cls_tbl$lipid_class[i], cls_tbl$n[i]))

  baseline_cond    <- paste0(BASELINE_TEMP, "_", BASELINE_TIME)
  baseline_samples <- meta %>% filter(condition == baseline_cond) %>% nrow()
  if (baseline_samples != N_REPLICATES)
    stop(sprintf("Expected %d baseline samples, found %d.", N_REPLICATES, baseline_samples))
  message("Baseline: ", baseline_cond, " - ", baseline_samples, " samples (OK)")

  rep_check <- meta %>%
    filter(condition != baseline_cond) %>%
    count(condition) %>%
    filter(n != N_REPLICATES)
  if (nrow(rep_check) > 0)
    warning("Some non-baseline conditions do not have exactly ", N_REPLICATES,
            " replicates:\n",
            paste(rep_check$condition, rep_check$n, sep = "=", collapse = ", "))
  else
    message("All non-baseline conditions have n=", N_REPLICATES, " replicates (OK)")

  list(long_df = long_df, meta = meta, lipid_order = lipid_order,
       lipid_class_map = lipid_class_map,
       lipid_display_labels = lipid_display_labels)
}

# =============================================================================
# 3. FOLD CHANGE
# =============================================================================

# log2FC = log2(median(condition) / median(baseline)) on raw peak areas.
compute_log2fc_matrix <- function(long_df) {
  baseline_cond <- paste0(BASELINE_TEMP, "_", BASELINE_TIME)

  medians <- long_df %>%
    group_by(lipid, condition) %>%
    summarise(median_val = median(raw_value), .groups = "drop")

  baseline_med <- medians %>%
    filter(condition == baseline_cond) %>%
    select(lipid, baseline_median = median_val)

  fc_wide <- medians %>%
    filter(condition != baseline_cond) %>%
    left_join(baseline_med, by = "lipid") %>%
    mutate(log2fc = log2(median_val / baseline_median)) %>%
    select(lipid, condition, log2fc) %>%
    pivot_wider(names_from = condition, values_from = log2fc) %>%
    arrange(lipid)

  message("log2FC matrix: ", nrow(fc_wide), " features x ",
          ncol(fc_wide) - 1, " conditions (baseline excluded)")
  fc_wide
}

order_fc_columns <- function(fc_wide) {
  panel_conds   <- unlist(HEATMAP_PANELS, use.names = FALSE)
  ordered_conds <- intersect(panel_conds, names(fc_wide))
  fc_wide %>% select(lipid, all_of(ordered_conds))
}

# =============================================================================
# 4. STATISTICS - per feature, per comparison group
# =============================================================================

run_kruskal_wallis_groups <- function(long_df, comparison_groups) {
  baseline_cond <- paste0(BASELINE_TEMP, "_", BASELINE_TIME)
  message("Running per-feature Kruskal-Wallis tests (one test per group)...")

  kw_list <- lapply(names(comparison_groups), function(grp) {
    grp_conds <- c(baseline_cond, comparison_groups[[grp]])
    sub_df    <- long_df %>% filter(condition %in% grp_conds)

    result <- sub_df %>%
      group_by(lipid) %>%
      summarise(
        kw_p = tryCatch(kruskal.test(log10_value ~ condition)$p.value,
                        error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      mutate(
        kw_p_adj = p.adjust(kw_p, method = "BH"),
        sig_raw  = kw_p     < 0.05,
        sig_BH   = kw_p_adj < SIG_THRESHOLD
      )

    message(sprintf("  [%s] %d / %d features significant at BH<%.2f",
                    grp, sum(result$sig_BH, na.rm = TRUE), nrow(result), SIG_THRESHOLD))

    names(result)[-1] <- paste0(names(result)[-1], "_", grp)
    result
  })

  kw_wide  <- Reduce(function(a, b) left_join(a, b, by = "lipid"), kw_list)
  sig_cols <- paste0("sig_BH_", names(comparison_groups))
  kw_wide %>%
    mutate(sig_any_BH = rowSums(across(all_of(sig_cols)), na.rm = TRUE) > 0)
}

# FSA::dunnTest() prints its full result table through message conditions on
# every call; suppressMessages() keeps long runs fast and logs readable.
run_dunn_posthoc_groups <- function(long_df, comparison_groups) {
  baseline_cond <- paste0(BASELINE_TEMP, "_", BASELINE_TIME)
  message("Running per-feature Dunn's post-hoc tests by group (may take several minutes)...")

  all_res <- lapply(names(comparison_groups), function(grp) {
    message(sprintf("  Group: %s", grp))
    grp_conds <- c(baseline_cond, comparison_groups[[grp]])
    sub_df    <- long_df %>% filter(condition %in% grp_conds)
    lipids    <- intersect(levels(long_df$lipid), unique(as.character(sub_df$lipid)))
    pb        <- txtProgressBar(min = 0, max = length(lipids), style = 3)

    grp_res <- vector("list", length(lipids))
    for (i in seq_along(lipids)) {
      lip    <- lipids[i]
      lip_df <- sub_df %>% filter(lipid == lip)
      tryCatch({
        dt <- suppressMessages(
          dunnTest(log10_value ~ condition, data = lip_df, method = "bh")$res
        )
        grp_res[[i]] <- dt %>%
          filter(str_detect(Comparison, fixed(baseline_cond))) %>%
          mutate(
            lipid        = lip,
            group        = grp,
            sig_star     = case_when(P.adj   < 0.001 ~ "***", P.adj   < 0.01 ~ "**",
                                     P.adj   < 0.05  ~ "*",   TRUE            ~ ""),
            sig_star_raw = case_when(P.unadj < 0.001 ~ "***", P.unadj < 0.01 ~ "**",
                                     P.unadj < 0.05  ~ "*",   TRUE            ~ "")
          )
      }, error = function(e) { grp_res[[i]] <<- NULL })
      setTxtProgressBar(pb, i)
    }
    close(pb)
    bind_rows(grp_res)
  })

  bind_rows(all_res) %>%
    select(lipid, group, Comparison, Z, P.unadj, P.adj, sig_star, sig_star_raw)
}

# =============================================================================
# 5. STATISTICS - per lipid class, per comparison group
# =============================================================================
# The class tests pool every replicate value of every lipid in the class.
# `sig_lipids` restricts the pool to a subset (significant or top-N lipids).

run_kruskal_wallis_class_groups <- function(long_df, comparison_groups, sig_lipids = NULL) {
  df            <- if (!is.null(sig_lipids)) long_df %>% filter(lipid %in% sig_lipids) else long_df
  baseline_cond <- paste0(BASELINE_TEMP, "_", BASELINE_TIME)
  message("Running per-class Kruskal-Wallis tests by group...")

  kw_list <- lapply(names(comparison_groups), function(grp) {
    grp_conds <- c(baseline_cond, comparison_groups[[grp]])
    sub_df    <- df %>% filter(condition %in% grp_conds)

    result <- sub_df %>%
      group_by(lipid_class) %>%
      summarise(
        kw_p = tryCatch(kruskal.test(log10_value ~ condition)$p.value,
                        error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      mutate(
        kw_p_adj = p.adjust(kw_p, method = "BH"),
        sig_raw  = kw_p     < 0.05,
        sig_BH   = kw_p_adj < SIG_THRESHOLD
      )

    names(result)[-1] <- paste0(names(result)[-1], "_", grp)
    result
  })

  kw_wide  <- Reduce(function(a, b) left_join(a, b, by = "lipid_class"), kw_list)
  sig_cols <- paste0("sig_BH_", names(comparison_groups))
  kw_wide %>%
    mutate(sig_any_BH = rowSums(across(all_of(sig_cols)), na.rm = TRUE) > 0)
}

run_dunn_posthoc_class_groups <- function(long_df, comparison_groups, sig_lipids = NULL) {
  df            <- if (!is.null(sig_lipids)) long_df %>% filter(lipid %in% sig_lipids) else long_df
  baseline_cond <- paste0(BASELINE_TEMP, "_", BASELINE_TIME)
  message("Running per-class Dunn's post-hoc tests by group...")

  all_res <- lapply(names(comparison_groups), function(grp) {
    grp_conds <- c(baseline_cond, comparison_groups[[grp]])
    sub_df    <- df %>% filter(condition %in% grp_conds)
    classes   <- unique(as.character(sub_df$lipid_class))
    classes   <- classes[!is.na(classes)]

    lapply(classes, function(cls) {
      cls_df <- sub_df %>% filter(lipid_class == cls)
      if (n_distinct(cls_df$condition) < 2) return(NULL)
      tryCatch({
        dt <- suppressMessages(
          dunnTest(log10_value ~ condition, data = cls_df, method = "bh")$res
        )
        dt %>%
          filter(str_detect(Comparison, fixed(baseline_cond))) %>%
          mutate(
            lipid_class  = cls,
            group        = grp,
            sig_star     = case_when(P.adj   < 0.001 ~ "***", P.adj   < 0.01 ~ "**",
                                     P.adj   < 0.05  ~ "*",   TRUE            ~ ""),
            sig_star_raw = case_when(P.unadj < 0.001 ~ "***", P.unadj < 0.01 ~ "**",
                                     P.unadj < 0.05  ~ "*",   TRUE            ~ "")
          )
      }, error = function(e) NULL)
    }) %>% bind_rows()
  }) %>% bind_rows()

  all_res %>%
    select(lipid_class, group, Comparison, Z, P.unadj, P.adj, sig_star, sig_star_raw)
}

# Within each class keep the `n` lipids with the smallest BH-adjusted KW p-value
# of TOP_N_RANK_GROUP (that group only; no combined score). Classes with fewer
# than n lipids keep all of them; lipids whose test failed (NA) rank last.
select_top_lipids_per_class <- function(kw_results, lipid_class_map, n = 10) {
  rank_col <- paste0("kw_p_adj_", TOP_N_RANK_GROUP)

  if (!rank_col %in% names(kw_results)) {
    warning("Ranking column '", rank_col, "' not found in the KW results - ",
            "is '", TOP_N_RANK_GROUP, "' a name in COMPARISON_GROUPS? ",
            "Skipping top-N-per-class selection.")
    return(character(0))
  }

  sel <- kw_results %>%
    transmute(lipid  = as.character(lipid),
              rank_p = .data[[rank_col]]) %>%
    left_join(lipid_class_map %>% mutate(lipid = as.character(lipid)),
              by = "lipid") %>%
    filter(!is.na(lipid_class)) %>%
    group_by(lipid_class) %>%
    arrange(rank_p, .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup()

  message(sprintf("Top-%d-per-class selection ranked on %s: %d lipids kept across %d classes",
                  n, rank_col, nrow(sel), n_distinct(sel$lipid_class)))
  sel$lipid
}

# =============================================================================
# 6. HEATMAPS
# =============================================================================

# Lookup env: key "<feature or class>|<condition>" -> asterisk string.
build_ast_env <- function(dunn_df, sig_col, key_col) {
  baseline_cond <- paste0(BASELINE_TEMP, "_", BASELINE_TIME)
  ast_map <- dunn_df %>%
    filter(.data[[sig_col]] != "") %>%
    mutate(
      cond_in_comp = case_when(
        str_starts(Comparison, paste0(baseline_cond, " - ")) ~
          str_replace(Comparison, paste0("^", baseline_cond, " - "), ""),
        str_ends(Comparison, paste0(" - ", baseline_cond)) ~
          str_replace(Comparison, paste0(" - ", baseline_cond, "$"), ""),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(cond_in_comp)) %>%
    mutate(key  = paste0(as.character(.data[[key_col]]), "|", cond_in_comp),
           star = .data[[sig_col]])

  env <- new.env(hash = TRUE, parent = emptyenv())
  for (k in seq_len(nrow(ast_map))) env[[ast_map$key[k]]] <- ast_map$star[k]
  env
}

panel_display_title <- function(panel_name) {
  if (panel_name %in% names(PANEL_DISPLAY_LABELS)) PANEL_DISPLAY_LABELS[[panel_name]] else panel_name
}

make_cell_fn <- function(ast_env, panel_conds, row_keys, fontsize) {
  if (is.null(ast_env)) return(NULL)
  force(panel_conds); force(row_keys); force(fontsize)
  function(j, i, x, y, width, height, fill) {
    ast <- ast_env[[paste0(row_keys[i], "|", panel_conds[j])]]
    if (!is.null(ast) && nchar(ast) > 0)
      grid.text(ast, x, y, gp = gpar(fontsize = fontsize, col = "black", fontface = "bold"))
  }
}

# Per-feature heatmap: rows grouped by lipid class in input order, no clustering.
make_lipid_heatmap <- function(fc_matrix, lipid_class_vec, row_display_labels,
                               dunn_df = NULL, sig_col = "sig_star") {

  mat <- fc_matrix %>% select(-lipid) %>% as.matrix()
  rownames(mat) <- as.character(fc_matrix$lipid)
  mat_clipped   <- pmin(pmax(mat, HEATMAP_LIMITS[1]), HEATMAP_LIMITS[2])
  col_fun       <- colorRamp2(c(HEATMAP_LIMITS[1], 0, HEATMAP_LIMITS[2]),
                              c("blue", "white", "red"))
  all_conds     <- colnames(mat)

  class_levels <- unique(as.character(lipid_class_vec))
  row_split    <- factor(as.character(lipid_class_vec), levels = class_levels)
  left_ann     <- rowAnnotation(
    Class = anno_block(gp = gpar(fill = "grey80", col = "grey50"),
                       labels = class_levels,
                       labels_gp = gpar(fontsize = 8, fontface = "bold"))
  )

  ast_env <- if (!is.null(dunn_df)) build_ast_env(dunn_df, sig_col, "lipid") else NULL

  valid_panels <- which(sapply(HEATMAP_PANELS,
                               function(pc) length(intersect(pc, all_conds)) > 0))
  last_pi <- max(valid_panels)
  hm_list <- NULL

  for (pi in seq_along(HEATMAP_PANELS)) {
    panel_name  <- names(HEATMAP_PANELS)[pi]
    panel_conds <- intersect(HEATMAP_PANELS[[pi]], all_conds)
    if (length(panel_conds) == 0) next

    sub_mat <- mat_clipped[, panel_conds, drop = FALSE]
    show_rn <- (pi == last_pi)

    ht <- Heatmap(
      sub_mat,
      name                  = paste0("log2FC_", pi),
      heatmap_legend_param  = list(title = "log2FC"),
      col                   = col_fun,
      cluster_rows          = FALSE,
      cluster_columns       = FALSE,
      row_split             = row_split,
      cluster_row_slices    = FALSE,
      row_title_gp          = gpar(fontsize = 0),
      column_title          = panel_display_title(panel_name),
      column_title_gp       = gpar(fontsize = 11, fontface = "bold"),
      column_labels         = str_replace(panel_conds, "^[^_]+_", ""),
      column_names_gp       = gpar(fontsize = 8),
      column_names_rot      = 0,
      column_names_centered = TRUE,
      show_row_names        = show_rn,
      row_labels            = if (show_rn) row_display_labels else character(0),
      row_names_gp          = gpar(fontsize = 7),
      row_names_side        = "right",
      left_annotation       = if (pi == 1) left_ann else NULL,
      show_heatmap_legend   = (pi == last_pi),
      cell_fun              = make_cell_fn(ast_env, panel_conds, rownames(sub_mat), 6),
      border                = TRUE,
      use_raster            = FALSE
    )
    hm_list <- if (is.null(hm_list)) ht else hm_list + ht
  }
  hm_list
}

# Class-summary heatmap: one row per lipid class (median log2FC of its lipids).
make_class_summary_heatmap <- function(class_fc_matrix, row_labels,
                                       dunn_class_df = NULL, sig_col = "sig_star") {
  mat <- class_fc_matrix %>%
    mutate(lipid_class = as.character(lipid_class)) %>%
    column_to_rownames("lipid_class") %>%
    as.matrix()
  mat_clipped <- pmin(pmax(mat, HEATMAP_LIMITS[1]), HEATMAP_LIMITS[2])
  col_fun     <- colorRamp2(c(HEATMAP_LIMITS[1], 0, HEATMAP_LIMITS[2]),
                            c("blue", "white", "red"))
  all_conds   <- colnames(mat)
  row_labels  <- unname(as.character(row_labels))

  ast_env <- if (!is.null(dunn_class_df)) build_ast_env(dunn_class_df, sig_col, "lipid_class") else NULL

  valid_panels <- which(sapply(HEATMAP_PANELS,
                               function(pc) length(intersect(pc, all_conds)) > 0))
  last_pi <- max(valid_panels)
  hm_list <- NULL

  for (pi in seq_along(HEATMAP_PANELS)) {
    panel_name  <- names(HEATMAP_PANELS)[pi]
    panel_conds <- intersect(HEATMAP_PANELS[[pi]], all_conds)
    if (length(panel_conds) == 0) next

    sub_mat <- mat_clipped[, panel_conds, drop = FALSE]

    ht <- Heatmap(
      sub_mat,
      name                  = paste0("log2FC_", pi),
      heatmap_legend_param  = list(title = "log2FC"),
      col                   = col_fun,
      cluster_rows          = FALSE,
      cluster_columns       = FALSE,
      column_title          = panel_display_title(panel_name),
      column_title_gp       = gpar(fontsize = 11, fontface = "bold"),
      column_labels         = str_replace(panel_conds, "^[^_]+_", ""),
      column_names_gp       = gpar(fontsize = 9),
      column_names_rot      = 0,
      column_names_centered = TRUE,
      show_row_names        = (pi == 1),
      row_labels            = if (pi == 1) row_labels else character(0),
      row_names_gp          = gpar(fontsize = 9),
      row_names_side        = "left",
      show_heatmap_legend   = (pi == last_pi),
      cell_fun              = make_cell_fn(ast_env, panel_conds, rownames(sub_mat), 8),
      border                = TRUE,
      use_raster            = FALSE
    )
    hm_list <- if (is.null(hm_list)) ht else hm_list + ht
  }
  hm_list
}

# =============================================================================
# 7. OUTPUT HELPERS
# =============================================================================

safe_write_table <- function(df, path_xlsx) {
  tryCatch({
    write.xlsx(df, path_xlsx, overwrite = TRUE)
    message("Wrote: ", path_xlsx)
  }, error = function(e) {
    path_csv <- str_replace(path_xlsx, "\\.xlsx$", ".csv")
    message("XLSX write failed (", conditionMessage(e),
            "); falling back to CSV: ", path_csv)
    write.csv(df, path_csv, row.names = FALSE)
  })
}

write_heatmap_pdf <- function(hm, path, width, height) {
  pdf(path, width = width, height = height)
  draw(hm, heatmap_legend_side = "right", merge_legend = TRUE)
  dev.off()
  message("Wrote: ", basename(path))
}

# =============================================================================
# 8. MAIN PIPELINE
# =============================================================================

main <- function() {
  tbl_dir <- file.path(OUTPUT_DIR, "tables")
  fig_dir <- file.path(OUTPUT_DIR, "figures")
  dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  topn <- paste0("top", TOP_N_PER_CLASS)

  # ---- Load -------------------------------------------------------------------
  loaded               <- load_data(INPUT_FILE, exclude_oxidized = EXCLUDE_OXIDIZED)
  long_df              <- loaded$long_df
  meta                 <- loaded$meta
  lipid_order          <- loaded$lipid_order
  lipid_class_map      <- loaded$lipid_class_map
  lipid_display_labels <- loaded$lipid_display_labels

  safe_write_table(meta, file.path(tbl_dir, "sample_metadata.xlsx"))

  log10_wide <- long_df %>%
    select(sample_id, lipid, log10_value) %>%
    pivot_wider(names_from = lipid, values_from = log10_value)
  safe_write_table(log10_wide, file.path(tbl_dir, "log10_transformed_data.xlsx"))

  # ---- Fold change ----------------------------------------------------------
  fc_wide <- order_fc_columns(compute_log2fc_matrix(long_df))
  safe_write_table(as.data.frame(fc_wide),
                   file.path(tbl_dir, "log2FC_matrix_per_feature.xlsx"))

  # ---- Per-feature statistics -----------------------------------------------
  kw_results <- run_kruskal_wallis_groups(long_df, COMPARISON_GROUPS)
  dunn_df    <- run_dunn_posthoc_groups(long_df, COMPARISON_GROUPS)
  safe_write_table(kw_results, file.path(tbl_dir, "kruskal_wallis_per_feature.xlsx"))
  safe_write_table(dunn_df,    file.path(tbl_dir, "dunn_posthoc_per_feature.xlsx"))

  sig_features <- kw_results %>% filter(sig_any_BH) %>% pull(lipid)
  message(sprintf("%d / %d features significant in >=1 group (KW BH<%.2f)",
                  length(sig_features), nrow(kw_results), SIG_THRESHOLD))

  lipid_class_vec_all <- lipid_class_map %>%
    mutate(lipid_class = factor(lipid_class, levels = levels(long_df$lipid_class))) %>%
    pull(lipid_class)

  sig_idx             <- which(lipid_order %in% as.character(sig_features))
  lipid_class_vec_sig <- lipid_class_vec_all[sig_idx]
  display_all         <- lipid_display_labels
  display_sig         <- lipid_display_labels[sig_idx]
  fc_wide_sig         <- fc_wide %>% filter(lipid %in% sig_features)

  # ---- Class-level analysis -------------------------------------------------
  build_class_fc <- function(lipids_subset = NULL) {
    fc_src    <- if (is.null(lipids_subset)) fc_wide else fc_wide %>% filter(lipid %in% lipids_subset)
    cond_cols <- setdiff(names(fc_src), "lipid")
    fc_src %>%
      mutate(lipid_class = factor(extract_lipid_class(as.character(lipid)),
                                  levels = levels(long_df$lipid_class))) %>%
      group_by(lipid_class) %>%
      summarise(across(all_of(cond_cols), \(x) median(x, na.rm = TRUE)), .groups = "drop") %>%
      filter(!is.na(lipid_class)) %>%
      arrange(lipid_class)
  }

  class_fc_all <- build_class_fc()
  class_fc_sig <- build_class_fc(sig_features)
  safe_write_table(class_fc_all, file.path(tbl_dir, "log2FC_matrix_per_class.xlsx"))
  safe_write_table(class_fc_sig, file.path(tbl_dir, "log2FC_matrix_per_class_significant.xlsx"))

  kw_class_all   <- run_kruskal_wallis_class_groups(long_df, COMPARISON_GROUPS)
  dunn_class_all <- run_dunn_posthoc_class_groups(long_df, COMPARISON_GROUPS)
  kw_class_sig   <- run_kruskal_wallis_class_groups(long_df, COMPARISON_GROUPS, sig_features)
  dunn_class_sig <- run_dunn_posthoc_class_groups(long_df, COMPARISON_GROUPS, sig_features)
  safe_write_table(kw_class_all,   file.path(tbl_dir, "kruskal_wallis_per_class.xlsx"))
  safe_write_table(dunn_class_all, file.path(tbl_dir, "dunn_posthoc_per_class.xlsx"))

  make_row_labels <- function(class_fc_mat, counts_df) {
    classes <- as.character(class_fc_mat$lipid_class)
    sapply(classes, function(cls) {
      n <- counts_df %>% filter(as.character(lipid_class) == cls) %>% pull(n)
      if (length(n) == 0) cls else paste0(cls, " (n=", n, ")")
    }, USE.NAMES = FALSE)
  }

  counts_all <- long_df %>% distinct(lipid, lipid_class) %>% count(lipid_class)
  counts_sig <- long_df %>% filter(lipid %in% sig_features) %>%
    distinct(lipid, lipid_class) %>% count(lipid_class)
  row_labels_all <- make_row_labels(class_fc_all, counts_all)
  row_labels_sig <- make_row_labels(class_fc_sig, counts_sig)

  # ---- Top-N-per-class family -----------------------------------------------
  top_features   <- select_top_lipids_per_class(kw_results, lipid_class_map, n = TOP_N_PER_CLASS)
  class_fc_top   <- NULL
  if (length(top_features) > 0) {
    rank_col <- paste0("kw_p_adj_", TOP_N_RANK_GROUP)

    top_detail <- lipid_class_map %>%
      mutate(lipid = as.character(lipid)) %>%
      filter(lipid %in% top_features) %>%
      left_join(kw_results %>% mutate(lipid = as.character(lipid)) %>%
                  select(lipid, all_of(rank_col)),
                by = "lipid") %>%
      mutate(lipid_class = factor(lipid_class, levels = levels(long_df$lipid_class))) %>%
      arrange(lipid_class, .data[[rank_col]])
    safe_write_table(top_detail, file.path(tbl_dir, paste0(topn, "_lipids_per_class.xlsx")))

    class_fc_top <- build_class_fc(top_features)
    safe_write_table(class_fc_top,
                     file.path(tbl_dir, paste0("log2FC_matrix_per_class_", topn, ".xlsx")))

    kw_class_top   <- run_kruskal_wallis_class_groups(long_df, COMPARISON_GROUPS, top_features)
    dunn_class_top <- run_dunn_posthoc_class_groups(long_df, COMPARISON_GROUPS, top_features)
    safe_write_table(kw_class_top,
                     file.path(tbl_dir, paste0("kruskal_wallis_per_class_", topn, ".xlsx")))
    safe_write_table(dunn_class_top,
                     file.path(tbl_dir, paste0("dunn_posthoc_per_class_", topn, ".xlsx")))

    counts_top <- long_df %>% filter(lipid %in% top_features) %>%
      distinct(lipid, lipid_class) %>% count(lipid_class)
    row_labels_top <- make_row_labels(class_fc_top, counts_top)
  } else {
    message("Top-N-per-class selection returned no lipids - skipping the top-N class family.")
  }

  # ---- Figures --------------------------------------------------------------
  pdf_w     <- 18
  pdf_h_all <- max(8, nrow(fc_wide)     * 0.12)
  pdf_h_sig <- max(8, nrow(fc_wide_sig) * 0.12)
  pdf_h_cls <- 7

  message("\n-- Per-feature heatmaps --")
  lhm <- function(fc, cls_vec, disp, dunn = NULL, sc = "sig_star")
    make_lipid_heatmap(fc, cls_vec, disp, dunn_df = dunn, sig_col = sc)

  write_heatmap_pdf(lhm(fc_wide, lipid_class_vec_all, display_all),
    file.path(fig_dir, "heatmap_all_features.pdf"), pdf_w, pdf_h_all)
  write_heatmap_pdf(lhm(fc_wide, lipid_class_vec_all, display_all, dunn_df, "sig_star"),
    file.path(fig_dir, "heatmap_all_features_significance_BH.pdf"), pdf_w, pdf_h_all)
  write_heatmap_pdf(lhm(fc_wide, lipid_class_vec_all, display_all, dunn_df, "sig_star_raw"),
    file.path(fig_dir, "heatmap_all_features_significance_uncorrected.pdf"), pdf_w, pdf_h_all)

  write_heatmap_pdf(lhm(fc_wide_sig, lipid_class_vec_sig, display_sig),
    file.path(fig_dir, "heatmap_significant_features.pdf"), pdf_w, pdf_h_sig)
  write_heatmap_pdf(lhm(fc_wide_sig, lipid_class_vec_sig, display_sig, dunn_df, "sig_star"),
    file.path(fig_dir, "heatmap_significant_features_significance_BH.pdf"), pdf_w, pdf_h_sig)
  write_heatmap_pdf(lhm(fc_wide_sig, lipid_class_vec_sig, display_sig, dunn_df, "sig_star_raw"),
    file.path(fig_dir, "heatmap_significant_features_significance_uncorrected.pdf"), pdf_w, pdf_h_sig)

  message("\n-- Class-summary heatmaps --")
  chm <- function(fc_cls, rl, dunn_cls = NULL, sc = "sig_star")
    make_class_summary_heatmap(fc_cls, rl, dunn_class_df = dunn_cls, sig_col = sc)

  write_heatmap_pdf(chm(class_fc_all, row_labels_all),
    file.path(fig_dir, "heatmap_class_all_features.pdf"), pdf_w, pdf_h_cls)
  write_heatmap_pdf(chm(class_fc_all, row_labels_all, dunn_class_all, "sig_star"),
    file.path(fig_dir, "heatmap_class_all_features_significance_BH.pdf"), pdf_w, pdf_h_cls)
  write_heatmap_pdf(chm(class_fc_all, row_labels_all, dunn_class_all, "sig_star_raw"),
    file.path(fig_dir, "heatmap_class_all_features_significance_uncorrected.pdf"), pdf_w, pdf_h_cls)

  if (nrow(class_fc_sig) > 0) {
    write_heatmap_pdf(chm(class_fc_sig, row_labels_sig),
      file.path(fig_dir, "heatmap_class_significant_features.pdf"), pdf_w, pdf_h_cls)
    write_heatmap_pdf(chm(class_fc_sig, row_labels_sig, dunn_class_sig, "sig_star"),
      file.path(fig_dir, "heatmap_class_significant_features_significance_BH.pdf"), pdf_w, pdf_h_cls)
    write_heatmap_pdf(chm(class_fc_sig, row_labels_sig, dunn_class_sig, "sig_star_raw"),
      file.path(fig_dir, "heatmap_class_significant_features_significance_uncorrected.pdf"), pdf_w, pdf_h_cls)
  } else {
    message("No significant lipids - skipping significant-class heatmaps.")
  }

  if (!is.null(class_fc_top) && nrow(class_fc_top) > 0) {
    write_heatmap_pdf(chm(class_fc_top, row_labels_top),
      file.path(fig_dir, paste0("heatmap_class_", topn, "_features.pdf")), pdf_w, pdf_h_cls)
    write_heatmap_pdf(chm(class_fc_top, row_labels_top, dunn_class_top, "sig_star"),
      file.path(fig_dir, paste0("heatmap_class_", topn, "_features_significance_BH.pdf")), pdf_w, pdf_h_cls)
    write_heatmap_pdf(chm(class_fc_top, row_labels_top, dunn_class_top, "sig_star_raw"),
      file.path(fig_dir, paste0("heatmap_class_", topn, "_features_significance_uncorrected.pdf")), pdf_w, pdf_h_cls)
  } else {
    message("No lipids selected for the top-N-per-class family - skipping those heatmaps.")
  }

  message("\n=== Pipeline complete. All outputs in: ", OUTPUT_DIR, " ===")
}

main()

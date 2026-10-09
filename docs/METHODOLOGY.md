# Methodology

This note describes the statistics behind `R/lipid_stability_analysis.R` and, in particular, how the three families of class-summary heatmaps differ.

## Pipeline overview

1. **Load and parse** the samples-in-rows table. Sample names are parsed with a regex built from `BASELINE_*`, `TEMP_ORDER`, `TIME_ORDER` and `FT_ORDER`. Lipid classes are derived from lipid names (`extract_lipid_class()`), then relabelled with `CLASS_MERGE_MAP`. If `EXCLUDE_OXIDIZED` is `TRUE`, all lipids whose class starts with `ox` are removed here, before any computation, so they are absent from every table and figure.
2. **Fold change** per lipid and condition: `log2(median(raw, condition) / median(raw, baseline))`. Raw (not log-transformed) values are used; the median is robust to outliers.
3. **Per-lipid screen**: Kruskal-Wallis test on log10 values, once per entry of `COMPARISON_GROUPS` (the baseline is added to every group). P-values are BH-adjusted across all lipids within a group. A lipid is flagged `sig_any_BH` if its adjusted p < `SIG_THRESHOLD` in at least one group.
4. **Per-lipid post-hoc**: `FSA::dunnTest(log10_value ~ condition, method = "bh")` per lipid and group; only baseline-versus-condition comparisons are kept. Both BH-adjusted (`sig_star`) and uncorrected (`sig_star_raw`) asterisks are stored.
5. **Class level**: the same Kruskal-Wallis and Dunn logic applied to all replicate values of all lipids of a class pooled together; the class colour is the median log2FC of its lipids.

Why separate groups? Running Kruskal-Wallis per comparison group prevents freeze-thaw conditions from diluting a storage-time signal and vice versa.

## Asterisks

`dunnTest()` returns two p-values per baseline-versus-condition comparison:

- `P.adj`: BH-adjusted *within that one class's (or lipid's) Dunn call for one comparison group*, i.e. across its own set of baseline comparisons, not across all classes or all lipids;
- `P.unadj`: the raw, uncorrected p-value.

`_significance_BH` figures draw stars from `P.adj`, `_significance_uncorrected` figures from `P.unadj`. Thresholds are identical: `*` p < 0.05, `**` p < 0.01, `***` p < 0.001. Dunn's test is run for every lipid and class regardless of whether the corresponding Kruskal-Wallis test was significant; the class-level Kruskal-Wallis results are stored in `kruskal_wallis_per_class*.xlsx` as a separate QC table and do not drive the stars.

## The three class-summary families

All class-summary heatmaps show one row per lipid class. They differ in exactly one respect: **which lipids are allowed to contribute to a class**, both for the colour value and for the class-level Dunn test. Everything downstream of that filter is identical.

| Family | Lipids pooled into each class | Figures |
|--------|-------------------------------|---------|
| **all** | every lipid in the class | `heatmap_class_all_features*` |
| **significant** | only lipids with `sig_any_BH == TRUE` from the per-lipid screen | `heatmap_class_significant_features*` |
| **top-N** | the `TOP_N_PER_CLASS` lipids per class with the smallest BH-adjusted Kruskal-Wallis p-value of `TOP_N_RANK_GROUP` | `heatmap_class_top10_features*` |

For each family and each class, the colour is the median log2FC over the included lipids, and the class-level Dunn test is run on the pooled raw replicate values of the included lipids. A class can therefore look different between families purely because different lipids were pooled, not because a different test was applied. For example, a class in which only 2 of 50 lipids passed the per-lipid screen is summarised by just those 2 lipids in the *significant* family.

### Top-N family

- Ranking uses **one** comparison group only (`TOP_N_RANK_GROUP`, default `Timepoint`). It deliberately does not use a combined or minimum-across-groups score: a group in which nothing responds would only add noise to the ordering.
- Lipids are ranked ascending, so lipids whose Kruskal-Wallis test failed (`NA`) are picked last and only if the class is small.
- Classes with fewer than `TOP_N_PER_CLASS` lipids keep all of them; there is no padding and no error. A class that ends up empty is skipped.
- The class-level Kruskal-Wallis and Dunn tests on the selected lipids still run for **all** comparison groups.
- `top10_lipids_per_class.xlsx` records exactly which lipids were selected, for QC. The file and figure names carry the actual value of `TOP_N_PER_CLASS` (`top5_...`, `top10_...`, ...).

## Design decisions

- **Median, not mean**, for fold change: robust to extreme values without a log transformation at that step.
- **log10** for the tests: a monotone transformation, so rank-based tests are unaffected by it, but it reduces skewness for the pooled class-level data.
- **Input order is locked**: lipids are factor-levelled in the order of the input columns; rows are never sorted or clustered.
- **`safe_write_table()`** falls back to CSV when an XLSX write fails (e.g. very wide tables).

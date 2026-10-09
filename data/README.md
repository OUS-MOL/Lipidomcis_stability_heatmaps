# data

Only **fictional** example data are tracked in this repository (`data/example/`). Real measurements are deliberately kept out: `.gitignore` ignores everything in `data/` except this file and `data/example/`. Put your own input files elsewhere (or in a git-ignored sub-folder such as `data/private/`).

## `data/example/`

| File | Content |
|------|---------|
| `example_POS.xlsx` | Synthetic Compound Discoverer-style export, positive mode: 105 compounds × 69 injections (60 design injections + blank, pooled-QC and bare-QC injections) |
| `example_NEG.xlsx` | Same for negative mode: 33 compounds × 69 injections |
| `example_combined_POS_NEG.xlsx` | Output of `R/combine_compound_discoverer.R`: 60 samples × 138 lipid features, **the input of the analysis script** |

All values are random numbers created by `scripts/generate_example_data.R` (fixed seed). Lipid names follow ordinary shorthand nomenclature, but the species, `m/z`, retention times and peak areas are invented, and the stability behaviour is artificial.

### Raw export layout (`example_POS.xlsx`, `example_NEG.xlsx`)

One compound per row:

| Column | Content |
|--------|---------|
| `Name` | Lipid annotation, e.g. `PC(16:0_18:1)` |
| `m/z`, `RT [min]` | Invented values |
| `Norm. Area: <injection>_fullMS_POS.raw (Fxx)` | One column per injection. Design injections are named `EXAMPLE_<temp>_<time>_<replicate>`, e.g. `EXAMPLE_-80_0h_1`, `EXAMPLE_RT_48h_3`, `EXAMPLE_-80_3xFT_2` |

The files also contain `Blank_processing_1`, `PQC_*` and bare `EXAMPLE_1…3` injections, which the combine script removes, and a few repeated lipid names (isomer-like duplicates) and PC species annotated in both polarities, so the de-duplication and `[POS]`/`[NEG]` tagging can be seen at work.

### Combined matrix layout (`example_combined_POS_NEG.xlsx`)

One sample per row: column `Sample` (e.g. `EXAMPLE_-80_0h_1`) followed by one column per lipid, with the ionisation mode appended to the name, e.g. `TG(O-50:1) [POS]`, `PE(20:4_22:6) [NEG]`. Repeated names within one mode are numbered before the tag, e.g. `TG(18:1_20:3_22:5) (1) [POS]` and `TG(18:1_20:3_22:5) (2) [POS]`.

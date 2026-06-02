# iDIFr <img src="man/figures/logo.png" align="right" height="139" alt="" />

**Intersectional Differential Item Functioning Analysis in R**

<!-- badges: start -->
[![R-CMD-check](https://github.com/thmsrgrs/iDIFr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thmsrgrs/iDIFr/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/iDIFr)](https://CRAN.R-project.org/package=iDIFr)
<!-- badges: end -->

`iDIFr` is an R package for detecting Differential Item Functioning (DIF)
using Logistic Regression, IRT Likelihood Ratio Tests, and model-based
recursive partitioning (MOB) — with first-class support for **intersectional
group designs** and built-in **Intersectional Contrast Analysis (ICA)**.

## Why iDIFr?

Most DIF packages focus on two-group comparisons along a single demographic
dimension. `iDIFr` is built around the idea that test-takers belong to multiple
groups simultaneously, and that DIF sometimes only appears at the *intersection*
of those identities.

Key features:

- **Intersectional group support** — define groups using `~ gender * nationality * age_band`
- **Effect sizes as first-class outputs** — results lead with Nagelkerke ΔR² and standardised chi, not just p-values
- **Three methods in one interface** — LR, LRT, and MOB with consistent output
- **Built-in ICA** — `ica = TRUE` classifies each item as *amplified*, *pure intersection*, *obscured*, or *none* by comparing single-variable and intersectional analyses
- **Transparent cell-size guidance** — `check_groups()` and `merge_groups()` help you manage sparse intersectional cells
- **Tidy output** — `tidy()` returns a flat data frame for use with `dplyr` and `ggplot2`

## Installation

```r
# From CRAN
install.packages("iDIFr")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("thmsrgrs/iDIFr")
```

## Quick start

```r
library(iDIFr)

# 1. Check your group structure first
check_groups(my_data, group = ~ gender * nationality * age_band)

# 2. Run DIF analysis — method selection is required
result <- idifr(
  data   = my_data,
  items  = 1:20,
  group  = ~ gender * nationality * age_band,
  method = c("LR", "LRT")
)

# 3. Explore results
print(result)                       # Flagged items with effect sizes
summary(result)                     # Full breakdown by method + concordance
plot(result)                        # Effect size heatmap
plot(result, type = "concordance")  # Method agreement
tidy(result)                        # Flat data frame
tidy(result, table = "direction")   # Group-level direction table
```

## Methods

| Argument | Method | Effect size | Best for |
|----------|--------|-------------|----------|
| `"LR"` | Logistic Regression | Nagelkerke ΔR² | General use, no IRT assumptions |
| `"LRT"` | IRT Likelihood Ratio Test | Standardised chi (df-scaled) | IRT-based programmes |
| `"MOB"` | Model-based recursive partitioning | Standardised score difference | Intersectional designs, exploratory |

## Intersectional Contrast Analysis (ICA)

Pass `ica = TRUE` to `idifr()` to run ICA automatically. After the main
analysis, `iDIFr` runs one additional `idifr()` per demographic variable and
classifies each item by comparing where it was flagged:

| Classification | Meaning |
|----------------|---------|
| `amplified` | Flagged in single-variable *and* intersectional runs |
| `pure_intersection` | Flagged only in the intersectional run |
| `obscured` | Flagged in a single-variable run but not intersectionally |
| `none` | Not flagged anywhere |

```r
result <- idifr(
  data   = my_data,
  items  = 1:20,
  group  = ~ gender * nationality * age_band,
  method = "LR",
  ica    = TRUE
)

print(result)                  # ICA section printed automatically
tidy(result, table = "ica")    # Flat ICA classification table
```

> **Note:** ICA runs N + 1 analyses without cross-analysis p-value correction.
> Interpret `pure_intersection` and `obscured` findings with caution in small
> samples.

## Effect size thresholds

`iDIFr` requires *both* statistical significance (after p-value adjustment)
*and* a meaningful effect size before flagging an item. This reduces false
positives in large samples.

| Method | Metric | Negligible | Moderate | Large |
|--------|--------|------------|----------|-------|
| LR (uniform) | Nagelkerke ΔR² | < .035 | .035–.070 | ≥ .070 |
| LR (non-uniform) | MAPPD | < .05 | .05–.10 | ≥ .10 |
| LRT (uniform) | Std. chi (df-scaled) | < 0.10×√(df/2) | 0.10–0.20×√(df/2) | ≥ 0.20×√(df/2) |
| LRT (non-uniform) | MAPPD | < .05 | .05–.10 | ≥ .10 |
| MOB | Std. score difference | < .35 | .35–.70 | ≥ .70 |

LRT thresholds are df-adjusted following Oshima et al. (1997) to maintain
equivalent sensitivity across designs with different numbers of groups.
The MOB threshold of 0.35 is intentionally conservative to avoid
over-detection in multigroup designs.

## Group management

```r
# Inspect cell sizes before analysis
check_groups(my_data, group = ~ gender * nationality * age_band)

# Merge sparse cells
grp <- check_groups(my_data, group = ~ gender * nationality * age_band)
merged_data <- merge_groups(
  grp,
  nationality = list("Other" = c("DE", "FR", "ES"))
)

# Merge multiple variables in one call
merged_data <- merge_groups(
  grp,
  nationality = list("Other" = c("DE", "FR")),
  age_band    = list("18-30" = c("18-24", "25-30"))
)

# Exclude groups below a minimum size at run time
result <- idifr(
  my_data, 1:20,
  group            = ~ gender * nationality * age_band,
  method           = "LR",
  exclude_below_min = TRUE,
  min_cell_size    = 50
)
```

## Simulating DIF data

`simulate_dif()` generates synthetic dichotomous item response data with
known DIF structure, including intersection-only DIF for validating `iDIFr`
on controlled data:

```r
# Standard DIF
dat <- simulate_dif(n_persons = 1000, n_items = 20, dif_items = c(3, 7))

# DIF confined to a single intersectional cell
dat_ix <- simulate_dif(
  n_persons     = 2000,
  n_items       = 20,
  dif_items     = c(5, 12),
  dif_effect    = 1.5,
  dif_structure = "intersection",
  dif_group     = list(group = "G1", nationality = "UK", age_band = "Young"),
  demo_vars     = list(nationality = c("UK", "DE", "FR"),
                       age_band    = c("Young", "Old")),
  seed          = 42
)

# Mixed DIF — some items standard, some intersectional
dat_mixed <- simulate_dif(
  n_persons     = 2000,
  n_items       = 20,
  dif_items     = list(standard = c(3, 7), intersection = c(12, 15)),
  dif_effect    = 1.0,
  dif_structure = "mixed",
  dif_group     = list(group = "G1", nationality = "UK", age_band = "Young"),
  demo_vars     = list(nationality = c("UK", "DE", "FR"),
                       age_band    = c("Young", "Old")),
  seed          = 42
)
```

## Citation

If you use `iDIFr` in published work, please cite:

```
Rogers, T. (2026). iDIFr: Intersectional Differential Item Functioning
Analysis in R. R package version 1.0.0.
https://CRAN.R-project.org/package=iDIFr
```

## Contributing

Bug reports and feature requests are welcome via
[GitHub Issues](https://github.com/thmsrgrs/iDIFr/issues).

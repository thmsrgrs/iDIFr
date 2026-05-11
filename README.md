# iDIFr <img src="man/figures/logo.png" align="right" height="139" alt="" />

**Intersectional Differential Item Functioning Analysis in R**

<!-- badges: start -->
[![R-CMD-check](https://github.com/username/iDIFr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/username/iDIFr/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/iDIFr)](https://CRAN.R-project.org/package=iDIFr)
<!-- badges: end -->

`iDIFr` is a user-friendly R package for detecting Differential Item
Functioning (DIF) using Logistic Regression, IRT Likelihood Ratio Tests, and
Random Forest structural-change tests — with first-class support for
**intersectional group designs**.

## Why iDIFr?

Most DIF packages focus on two-group comparisons along a single demographic
dimension. `iDIFr` is built around the idea that test-takers carry multiple
identities simultaneously, and that DIF sometimes only appears at the
*intersection* of those identities — for example, in the group "women from
non-English-speaking backgrounds aged 18–30", rather than in any single group.

Key features:

- **Intersectional group support** — define groups using `~ gender * nationality * age_band`
- **Effect sizes as first-class outputs** — results lead with Nagelkerke ΔR² and ΔCFI, not just p-values
- **Three methods in one interface** — LR, LRT, and RF with consistent output
- **Transparent cell-size guidance** — `check_groups()` and `merge_groups()` help you manage sparse intersectional cells
- **Tidy output** — `tidy()` returns a flat data frame for use with `dplyr` and `ggplot2`

## Installation

```r
# Development version from GitHub
# install.packages("remotes")
remotes::install_github("username/iDIFr")
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
print(result)          # Flagged items with effect sizes
summary(result)        # Full breakdown by method + concordance
plot(result)           # Effect size heatmap
plot(result, type = "concordance")  # Method agreement
tidy(result)           # Flat data frame
```

## Methods

| Argument | Method | Effect size | Best for |
|----------|--------|-------------|----------|
| `"LR"` | Logistic Regression | Nagelkerke ΔR² (ETS A/B/C) | General use, no IRT assumptions |
| `"LRT"` | IRT Likelihood Ratio Test | Std. chi / ΔRMSEA | IRT-based programmes |
| `"RF"` | Random Forest (structural change) | Std. score difference | Intersectional designs, no linearity assumption |

## Effect size thresholds

`iDIFr` requires *both* statistical significance (after p-value adjustment)
*and* a meaningful effect size before flagging an item. This reduces false
positives in large samples.

| Method | Negligible | Moderate | Large |
|--------|-----------|---------|-------|
| LR (ΔR²) | < .035 | .035–.070 | ≥ .070 |
| LRT (std. chi) | < 0.2 | 0.2–0.5 | ≥ .050 |
| RF (std. diff) | < .20 | .20–.50 | ≥ .50 |

## Citation

If you use `iDIFr` in published work, please cite:

```
Author (2026). iDIFr: Intersectional Differential Item Functioning Analysis
in R. R package version 0.1.0. https://github.com/username/iDIFr
```

## Contributing

Bug reports and feature requests are welcome via
[GitHub Issues](https://github.com/username/iDIFr/issues).

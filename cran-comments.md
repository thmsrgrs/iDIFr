## Resubmission

This is a resubmission addressing the CRAN reviewer's comments:

* Removed "in R" from package title
* Added single quotes around package and software names in DESCRIPTION
* Added method references to DESCRIPTION
* Added \value tags to plot.idifr(), print.idifr(), and summary.idifr()
* Replaced \dontrun{} with \donttest{} in examples
* Standardised LR effect size labels to Negligible/Moderate/Large

## Test environments
* Local Windows 10, R 4.4.0: 0 errors, 0 warnings, 1 note
* win-builder R-devel: 0 errors, 0 warnings, 1 note

## R CMD check results
0 errors | 0 warnings | 1 note

* NOTE: Possibly misspelled words in DESCRIPTION — these are author
  surnames (Thissen, Steinberg, Wainer, Swaminathan, Strobl, Kopf,
  Zeileis) and standard psychometric abbreviations (DIF, ICA).
  All are spelled correctly.
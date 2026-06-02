## Resubmission

This is a resubmission addressing the CRAN reviewer's comments:

* Removed "in R" from package title
* Added single quotes around package and software names in DESCRIPTION
* Added method references to DESCRIPTION
* Added \value tags to plot.idifr(), print.idifr(), and summary.idifr()
* Replaced \dontrun{} with \donttest{} in examples

## Test environments
* Local Windows 10, R 4.4.0: 0 errors, 0 warnings, 1 note
* win-builder R-devel: 0 errors, 0 warnings, 2 notes

## R CMD check results
0 errors | 0 warnings | 2 notes

* NOTE: "unable to verify current time" — unreachable time server,
  not a package issue.

* NOTE: DIF, ICA, and iDIFr flagged as possibly misspelled — these
  are standard psychometric abbreviations and the package name.
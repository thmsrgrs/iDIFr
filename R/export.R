# export.R — Export iDIFr results to Excel
#
# Writes an idifr result object to a formatted .xlsx workbook.
# Requires the 'openxlsx' package (listed in Suggests).

#' Export iDIFr results to Excel
#'
#' @description
#' Writes an `idifr` result object to a formatted `.xlsx` workbook.  Each
#' requested sheet is written as an Excel table so that column headers are
#' bold, filters are enabled, and values are properly typed.
#'
#' Only columns that are actually present in the result object are written;
#' columns listed in the per-method definitions that were not produced by the
#' current run are silently omitted.
#'
#' @param x An `idifr` object returned by [idifr()].
#' @param file Path to the output `.xlsx` file (character string).
#' @param sheets Character vector of sheet keys to include.  Valid keys:
#'   `"summary"`, `"lr"`, `"lrt"`, `"mob"`, `"direction"`, `"ica"`,
#'   `"groups"`.  Pass `NULL` (default) to include all available sheets.
#' @param overwrite Logical.  If `TRUE` (default) an existing file is
#'   silently overwritten.
#'
#' @return `x` invisibly (so the call can be piped).
#'
#' @examples
#' \dontrun{
#' result <- idifr(dat, 1:20, ~ gender * country, method = c("LR", "LRT"),
#'                 ica = TRUE)
#' export_results(result, "my_dif_results.xlsx")
#'
#' # Only summary and ICA sheets
#' export_results(result, "summary_only.xlsx", sheets = c("summary", "ica"))
#' }
#'
#' @export
export_results <- function(x, file,
                            sheets    = NULL,
                            overwrite = TRUE) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop(
      "Package 'openxlsx' is required for export_results(). ",
      "Install it with: install.packages('openxlsx')",
      call. = FALSE
    )
  }

  if (!inherits(x, "idifr")) {
    stop("`x` must be an idifr object returned by idifr().", call. = FALSE)
  }

  # ---- Column order definitions per method ----------------------------------
  # Only columns present in the actual result data are written (safe_select).
  # The order controls column display in Excel.

  lr_cols <- c(
    "item", "flagged", "dif_type",
    "chi_sq_uniform",    "df_uniform",    "p_uniform",
    "chi_sq_nonuniform", "df_nonuniform", "p_nonuniform",
    "chi_sq_interaction", "p_interaction",
    "p_overall", "p_adj", "p_adj_interaction",
    "delta_r2_uniform", "delta_r2_interaction",
    "delta_r2_omnibus", "delta_r2",
    "mappd", "nu_es", "nu_es_class", "ets_class"
  )

  lrt_cols <- c(
    "item", "flagged", "dif_type",
    "chi_sq", "df", "std_chi", "es_class",
    "chi_uniform",    "df_uniform",    "p_uniform_adj",
    "std_chi_uniform",    "es_class_uniform",
    "chi_nonuniform", "df_nonuniform", "p_nonuniform_adj",
    "std_chi_nonuniform", "mappd", "mappd_class",
    "p_overall", "p_adj"
  )

  mob_cols <- c(
    "item", "flagged", "dif_type",
    "dif_source", "split_variable", "split_depth", "split_vars_all",
    "std_diff", "es_class",
    "p_split", "p_overall", "p_adj"
  )

  safe_select <- function(df, cols) {
    df[, intersect(cols, names(df)), drop = FALSE]
  }

  all_results  <- x$results
  methods_run  <- unique(all_results$method)

  # ---- Per-method tabs -------------------------------------------------------

  tabs <- list()

  if ("LR" %in% methods_run) {
    tabs$lr <- safe_select(
      all_results[all_results$method == "LR", ],
      lr_cols
    )
  }

  if ("LRT" %in% methods_run) {
    tabs$lrt <- safe_select(
      all_results[all_results$method == "LRT", ],
      lrt_cols
    )
  }

  if ("MOB" %in% methods_run) {
    tabs$mob <- safe_select(
      all_results[all_results$method == "MOB", ],
      mob_cols
    )
  }

  # ---- Summary tab: key columns from every method side-by-side --------------

  sum_parts <- list()

  if ("LR" %in% methods_run) {
    lr_s <- safe_select(
      tabs$lr,
      c("item", "flagged", "dif_type", "delta_r2", "ets_class", "p_adj")
    )
    names(lr_s)[names(lr_s) != "item"] <-
      paste0("LR_", names(lr_s)[names(lr_s) != "item"])
    sum_parts$lr <- lr_s
  }

  if ("LRT" %in% methods_run) {
    lrt_s <- safe_select(
      tabs$lrt,
      c("item", "flagged", "dif_type", "std_chi", "es_class", "p_adj")
    )
    names(lrt_s)[names(lrt_s) != "item"] <-
      paste0("LRT_", names(lrt_s)[names(lrt_s) != "item"])
    sum_parts$lrt <- lrt_s
  }

  if ("MOB" %in% methods_run) {
    mob_s <- safe_select(
      tabs$mob,
      c("item", "flagged", "dif_type", "std_diff", "es_class", "p_adj")
    )
    names(mob_s)[names(mob_s) != "item"] <-
      paste0("MOB_", names(mob_s)[names(mob_s) != "item"])
    sum_parts$mob <- mob_s
  }

  if (!is.null(x$ica)) {
    # Use LR ICA if available, otherwise first method available
    ica_method <- if ("LR" %in% methods_run) "LR" else methods_run[1]
    ica_s <- x$ica[x$ica$method == ica_method,
                    c("item", "ica_class", "marginal_vars"),
                    drop = FALSE]
    names(ica_s)[names(ica_s) != "item"] <-
      paste0("ICA_", names(ica_s)[names(ica_s) != "item"])
    sum_parts$ica <- ica_s
  }

  summary_tab <- if (length(sum_parts) > 0) {
    Reduce(function(a, b) merge(a, b, by = "item", all = TRUE), sum_parts)
  } else {
    NULL
  }

  # ---- Determine which sheets are available ---------------------------------

  available <- "summary"
  if ("LR"  %in% methods_run)           available <- c(available, "lr")
  if ("LRT" %in% methods_run)           available <- c(available, "lrt")
  if ("MOB" %in% methods_run)           available <- c(available, "mob")
  if (!is.null(x$group_direction) &&
      nrow(x$group_direction) > 0)      available <- c(available, "direction")
  if (!is.null(x$ica))                  available <- c(available, "ica")
  if (!is.null(x$groups$cell_counts))   available <- c(available, "groups")

  if (is.null(sheets)) {
    sheets <- available
  } else {
    unknown <- setdiff(sheets, c("summary","lr","lrt","mob",
                                  "direction","ica","groups"))
    if (length(unknown) > 0) {
      warning(
        "Unknown sheet key(s) ignored: ", paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    sheets <- intersect(sheets, available)
  }

  if (length(sheets) == 0) {
    cli::cli_alert_warning(
      "No sheets to write — check that the requested sheets are available."
    )
    return(invisible(x))
  }

  # ---- Build workbook -------------------------------------------------------

  sheet_order <- c("summary", "lr", "lrt", "mob", "direction", "ica", "groups")

  wb <- openxlsx::createWorkbook()

  for (s in intersect(sheet_order, sheets)) {

    sheet_name <- switch(s,
      summary   = "Summary",
      lr        = "LR",
      lrt       = "LRT",
      mob       = "MOB",
      direction = "Direction",
      ica       = "ICA",
      groups    = "Groups"
    )

    dat <- switch(s,
      summary   = summary_tab,
      lr        = tabs$lr,
      lrt       = tabs$lrt,
      mob       = tabs$mob,
      direction = x$group_direction,
      ica       = x$ica,
      groups    = x$groups$cell_counts
    )

    if (!is.null(dat) && nrow(dat) > 0) {
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeDataTable(wb, sheet_name, dat,
                                tableName = paste0("tbl_", sheet_name))
    }
  }

  openxlsx::saveWorkbook(wb, file, overwrite = overwrite)

  cli::cli_alert_success("Results saved to {.path {file}}")

  invisible(x)
}

#' Run intersectional DIF analysis
#'
#' @description
#' The main entry point for `iDIFr`. Detects Differential Item Functioning (DIF)
#' using one or more statistical methods, with full support for intersectional
#' group structures defined by crossing multiple demographic variables.
#'
#' Effect sizes are reported alongside significance for all methods. Groups with
#' small cell sizes trigger a warning. Use `exclude_below_min` and
#' `fully_crossed` to control whether those groups are included in the analysis.
#'
#' @param data A data frame containing item responses and demographic variables.
#' @param items A numeric vector of column indices, or a character vector of
#'   column names, identifying the item response columns. Items must be
#'   dichotomously scored (0/1).
#' @param group A one-sided formula specifying the grouping variable(s).
#'   Use `~ var` for a single demographic (2+ groups) or `~ var1 * var2` for
#'   intersectional groups. Example: `~ gender * nationality * age_band`.
#' @param method A character vector specifying which DIF method(s) to use.
#'   Must be one or more of `"LR"` (Logistic Regression), `"LRT"` (IRT
#'   Likelihood Ratio Test), `"ID"` (Intersectional Decomposition), or `"RF"`
#'   (Random Forest structural-change test). No default -- the user must
#'   choose.
#' @param min_cell_size Minimum acceptable group size. Groups below this
#'   threshold trigger a warning. Also used as the crossing criterion when
#'   `exclude_below_min = TRUE` or `fully_crossed` is supplied. Default is 50.
#' @param exclude_below_min Logical. If `TRUE`, any intersectional group with
#'   fewer than `min_cell_size` respondents is excluded from the analysis
#'   entirely. If `FALSE` (default), all groups are included and small groups
#'   trigger a warning only.
#' @param fully_crossed A character vector of variable name(s). Only levels of
#'   the named variable(s) that are fully crossed -- meaning every intersectional
#'   cell for that level meets `min_cell_size` -- are included in the analysis.
#'   Respondents belonging to levels that are not fully crossed are excluded.
#'   Default is `NULL` (no crossing filter applied). Example:
#'   `fully_crossed = "nationality"` keeps only nationalities where every
#'   gender x age_band cell meets `min_cell_size`.
#' @param value_selection A named list for filtering specific values of
#'   demographic variables before analysis. Each element should be named after
#'   a grouping variable and contain a character vector of values to keep.
#'   Variables not mentioned are left unchanged (all values included). Default
#'   is `NULL`. Example: `value_selection = list(country = c("UK", "France"),
#'   age_band = c("Young", "Old"))`.
#' @param anchor A numeric or character vector identifying anchor items
#'   (items assumed to be DIF-free) for IRT scaling. If `NULL` (default),
#'   all items are used as anchors in the first pass.
#' @param alpha Significance level for DIF flagging. Default is `0.05`.
#' @param p_adjust Method for p-value adjustment across items. Passed to
#'   [stats::p.adjust()]. Default is `"BH"` (Benjamini-Hochberg). Use
#'   `"none"` to skip adjustment.
#' @param verbose Logical. If `TRUE` (default), prints progress and group
#'   information during the analysis.
#'
#' @return An object of class `idifr` containing:
#' \describe{
#'   \item{results}{A data frame with one row per item per method, including
#'     test statistics, p-values, adjusted p-values, effect sizes, and DIF
#'     classification (A/B/C for LR; negligible/moderate/large for LRT and ID).}
#'   \item{groups}{An `idifr_groups` object describing the group structure,
#'     cell sizes, and any small-cell warnings.}
#'   \item{method}{Character vector of methods used.}
#'   \item{call}{The matched call.}
#'   \item{items}{Character vector of item names analysed.}
#'   \item{alpha}{The significance level used.}
#'   \item{p_adjust}{The p-value adjustment method used.}
#'   \item{excluded_groups}{Character vector of group labels excluded by
#'     `exclude_below_min` or `fully_crossed`, or `NULL` if no exclusions.}
#'   \item{excluded_values}{Named list of value_selection filters applied,
#'     or `NULL` if none.}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic two-group analysis
#' result <- idifr(
#'   data   = my_data,
#'   items  = 1:20,
#'   group  = ~ gender,
#'   method = c("LR", "LRT")
#' )
#'
#' # Intersectional analysis, excluding groups below minimum
#' result <- idifr(
#'   data              = my_data,
#'   items             = 1:20,
#'   group             = ~ gender * nationality * age_band,
#'   method            = c("LR", "LRT"),
#'   min_cell_size     = 50,
#'   exclude_below_min = TRUE
#' )
#'
#' # Only include nationalities fully crossed across gender and age_band
#' result <- idifr(
#'   data          = my_data,
#'   items         = 1:20,
#'   group         = ~ gender * nationality * age_band,
#'   method        = "LR",
#'   fully_crossed = "nationality"
#' )
#'
#' # Restrict to specific variable values
#' result <- idifr(
#'   data            = my_data,
#'   items           = 1:20,
#'   group           = ~ gender * nationality * age_band,
#'   method          = "LR",
#'   value_selection = list(nationality = c("UK", "France"),
#'                          age_band    = c("Young", "Old"))
#' )
#' }
#'
#' @seealso [check_groups()] for exploring group structure before analysis;
#'   [group_details()] and [cross_details()] for full breakdowns;
#'   [merge_groups()] for combining sparse cells.
#'
#' @export
idifr <- function(data,
                  items,
                  group,
                  method,
                  min_cell_size     = 50,
                  exclude_below_min = FALSE,
                  fully_crossed     = NULL,
                  value_selection   = NULL,
                  anchor            = NULL,
                  alpha             = 0.05,
                  p_adjust          = "BH",
                  verbose           = TRUE) {

  call <- match.call()

  # --- Input validation -------------------------------------------------------

  .validate_inputs(data, items, group, method, alpha, p_adjust,
                   exclude_below_min, fully_crossed, value_selection)

  # --- Apply value_selection filter ------------------------------------------

  if (!is.null(value_selection)) {
    data <- .apply_value_selection(data, value_selection, group, verbose)
  }

  # --- Resolve item columns ---------------------------------------------------

  item_cols <- .resolve_items(data, items)

  # --- Construct groups -------------------------------------------------------

  groups <- .build_groups(data, group, min_cell_size, verbose)

  # --- Apply fully_crossed filter --------------------------------------------

  excluded_groups <- NULL

  if (!is.null(fully_crossed)) {
    filter_result  <- .apply_fully_crossed(data, groups, fully_crossed,
                                           min_cell_size, verbose)
    data           <- filter_result$data
    excluded_groups <- c(excluded_groups, filter_result$excluded)
    # Rebuild groups on filtered data
    groups <- .build_groups(data, group, min_cell_size, verbose = FALSE)
  }

  # --- Apply exclude_below_min filter ----------------------------------------

  if (isTRUE(exclude_below_min)) {
    filter_result   <- .apply_exclude_below_min(data, groups, min_cell_size,
                                                verbose)
    data            <- filter_result$data
    excluded_groups <- c(excluded_groups, filter_result$excluded)
    # Rebuild groups on filtered data
    groups <- .build_groups(data, group, min_cell_size, verbose = FALSE)
  }

  # --- Warn if small cells remain --------------------------------------------

  if (groups$has_small_cells && verbose) {
    cli::cli_alert_warning(
      "One or more groups have fewer than {min_cell_size} respondents."
    )
    cli::cli_alert_info(
      "Run {.code check_groups(data, group = {deparse(group)})} for a \\
       full breakdown. Use {.code exclude_below_min = TRUE} or \\
       {.code fully_crossed} to exclude small groups from the analysis."
    )
    cat("\n")
  }

  # --- Report exclusions -----------------------------------------------------

  if (!is.null(excluded_groups) && verbose) {
    n_excl <- length(unique(excluded_groups))
    cli::cli_alert_info(
      "{n_excl} group(s) excluded from analysis. \\
       See {.code result$excluded_groups} for details."
    )
    cat("\n")
  }

  # --- Detect design type -----------------------------------------------------

  design <- .detect_design(group, data)

  if (verbose) {
    cli::cli_alert_info(
      "Design: {design} ({groups$n_groups} group{?s})"
    )
    cat("\n")
  }

  # --- Dispatch to method modules --------------------------------------------

  results_list <- list()

  if ("LR" %in% method) {
    if (verbose) cli::cli_h2("Running Logistic Regression DIF (LR)")
    results_list[["LR"]] <- .run_lr(data, item_cols, groups, alpha,
                                    p_adjust, verbose)
  }

  if ("LRT" %in% method) {
    if (verbose) cli::cli_h2("Running IRT Likelihood Ratio Test (LRT)")
    results_list[["LRT"]] <- .run_lrt(data, item_cols, groups, anchor,
                                      alpha, p_adjust, verbose)
  }

  if ("ID" %in% method) {
    if (verbose) cli::cli_h2("Running Intersectional Decomposition (ID)")
    results_list[["ID"]] <- .run_id(data, item_cols, groups, alpha,
                                    p_adjust, verbose)
  }

  if ("RF" %in% method) {
    if (verbose) cli::cli_h2("Running Random Forest DIF (RF)")
    results_list[["RF"]] <- .run_rf(data, item_cols, groups, alpha,
                                    p_adjust, verbose)
  }

  # --- Combine results -------------------------------------------------------

  combined <- .combine_results(results_list, item_cols)

  # --- Build and return idifr object -----------------------------------------

  structure(
    list(
      results          = combined$results,
      group_direction  = combined$group_direction,
      groups           = groups,
      design           = design,
      method           = method,
      call             = call,
      items            = item_cols,
      alpha            = alpha,
      p_adjust         = p_adjust,
      excluded_groups  = if (length(excluded_groups) > 0) unique(excluded_groups) else NULL,
      excluded_values  = value_selection
    ),
    class = "idifr"
  )
}


# --- Internal: detect design type ---------------------------------------------

.detect_design <- function(group_formula, data) {
  vars     <- all.vars(group_formula)
  n_vars   <- length(vars)
  n_levels <- sapply(vars, function(v) length(unique(data[[v]])))
  total_groups <- prod(n_levels)

  if      (n_vars == 1 && total_groups == 2) "dichotomous"
  else if (n_vars == 1 && total_groups >  2) "multigroup"
  else                                        "intersectional"
}


# --- Internal validation ------------------------------------------------------

.validate_inputs <- function(data, items, group, method, alpha, p_adjust,
                              exclude_below_min, fully_crossed, value_selection) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (missing(method) || is.null(method)) {
    stop(
      "You must specify at least one method.\n",
      "  Choose from: \"LR\", \"LRT\", \"ID\", \"RF\"\n",
      "  Example: method = c(\"LR\", \"LRT\")",
      call. = FALSE
    )
  }

  valid_methods <- c("LR", "LRT", "ID", "RF")
  bad_methods <- setdiff(toupper(method), valid_methods)
  if (length(bad_methods) > 0) {
    stop(
      "Unknown method(s): ", paste(bad_methods, collapse = ", "), "\n",
      "  Valid options are: ", paste(valid_methods, collapse = ", "),
      call. = FALSE
    )
  }

  if (!inherits(group, "formula")) {
    stop(
      "`group` must be a formula. Example: group = ~ gender * nationality",
      call. = FALSE
    )
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a number between 0 and 1.", call. = FALSE)
  }

  valid_adjust <- c(stats::p.adjust.methods, "none")
  if (!p_adjust %in% valid_adjust) {
    stop(
      "`p_adjust` must be one of: ", paste(valid_adjust, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.logical(exclude_below_min) || length(exclude_below_min) != 1) {
    stop("`exclude_below_min` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(fully_crossed) && !is.character(fully_crossed)) {
    stop(
      "`fully_crossed` must be a character vector of variable name(s).\n",
      "  Example: fully_crossed = \"nationality\"",
      call. = FALSE
    )
  }

  if (!is.null(value_selection)) {
    if (!is.list(value_selection) || is.null(names(value_selection))) {
      stop(
        "`value_selection` must be a named list.\n",
        "  Example: value_selection = list(nationality = c(\"UK\", \"France\"))",
        call. = FALSE
      )
    }
    vars <- all.vars(group)
    bad_vs <- setdiff(names(value_selection), vars)
    if (length(bad_vs) > 0) {
      stop(
        "These `value_selection` variables are not in the group formula: ",
        paste(bad_vs, collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


# --- Internal: apply value_selection filter -----------------------------------

.apply_value_selection <- function(data, value_selection, group, verbose) {

  original_n <- nrow(data)

  for (var in names(value_selection)) {
    keep_vals <- value_selection[[var]]
    data <- data[as.character(data[[var]]) %in% as.character(keep_vals), ]
  }

  n_removed <- original_n - nrow(data)

  if (verbose && n_removed > 0) {
    filter_desc <- paste(
      mapply(function(v, vals) paste0(v, " = c(", paste(vals, collapse = ", "), ")"),
             names(value_selection), value_selection),
      collapse = "; "
    )
    cli::cli_alert_info(
      "value_selection: {n_removed} respondent(s) excluded. \\
       Filters applied: {filter_desc}"
    )
    cat("\n")
  }

  if (nrow(data) == 0) {
    stop(
      "`value_selection` filters removed all respondents. \\
       Check that the specified values exist in the data.",
      call. = FALSE
    )
  }

  data
}


# --- Internal: apply fully_crossed filter -------------------------------------

.apply_fully_crossed <- function(data, groups, fully_crossed, min_cell_size,
                                  verbose) {

  vars <- groups$vars
  tbl  <- groups$cell_table

  bad_vars <- setdiff(fully_crossed, vars)
  if (length(bad_vars) > 0) {
    stop(
      "These `fully_crossed` variables are not in the group formula: ",
      paste(bad_vars, collapse = ", "),
      call. = FALSE
    )
  }

  excluded_levels <- list()

  for (cv in fully_crossed) {

    cross_res <- .check_crossing(tbl, cv, vars, min_cell_size)
    not_crossed_levels <- cross_res$level[!cross_res$fully_crossed]

    if (length(not_crossed_levels) > 0) {
      excluded_levels[[cv]] <- as.character(not_crossed_levels)
      data <- data[!as.character(data[[cv]]) %in% as.character(not_crossed_levels), ]

      if (verbose) {
        cli::cli_alert_info(
          "fully_crossed: {length(not_crossed_levels)} level(s) of \\
           {.field {cv}} excluded (not fully crossed): \\
           {paste(as.character(not_crossed_levels), collapse = ', ')}"
        )
      }
    } else {
      if (verbose) {
        cli::cli_alert_success(
          "fully_crossed: all levels of {.field {cv}} are fully crossed. \\
           No exclusions."
        )
      }
    }
  }

  if (verbose && length(excluded_levels) > 0) cat("\n")

  # Build excluded group labels for reporting
  excluded_group_labels <- character(0)
  if (length(excluded_levels) > 0) {
    excl_rows <- Reduce(`|`, lapply(names(excluded_levels), function(cv) {
      as.character(groups$data[[cv]]) %in% excluded_levels[[cv]]
    }))
    excluded_group_labels <- unique(as.character(groups$group_vector[excl_rows]))
  }

  list(data = data, excluded = excluded_group_labels)
}


# --- Internal: apply exclude_below_min filter ---------------------------------

.apply_exclude_below_min <- function(data, groups, min_cell_size, verbose) {

  small <- groups$cell_counts[groups$cell_counts$n < min_cell_size, ".group_label"]

  if (length(small) == 0) {
    if (verbose) {
      cli::cli_alert_success(
        "exclude_below_min: no groups below minimum -- no exclusions made."
      )
      cat("\n")
    }
    return(list(data = data, excluded = character(0)))
  }

  data <- data[!as.character(groups$group_vector) %in% small, ]

  if (verbose) {
    cli::cli_alert_info(
      "exclude_below_min: {length(small)} group(s) with n < {min_cell_size} \\
       excluded: {paste(small, collapse = ', ')}"
    )
    cat("\n")
  }

  list(data = data, excluded = small)
}


# --- Resolve item columns -----------------------------------------------------

.resolve_items <- function(data, items) {
  if (is.numeric(items)) {
    if (any(items < 1) || any(items > ncol(data))) {
      stop("Some item indices are out of range for the supplied data frame.",
           call. = FALSE)
    }
    cols <- names(data)[items]
  } else if (is.character(items)) {
    missing_cols <- setdiff(items, names(data))
    if (length(missing_cols) > 0) {
      stop("These item columns were not found in `data`: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    cols <- items
  } else {
    stop("`items` must be a numeric vector of column indices or a character \\
          vector of column names.", call. = FALSE)
  }

  # Check items are dichotomous
  for (col in cols) {
    vals <- unique(stats::na.omit(data[[col]]))
    if (!all(vals %in% c(0, 1))) {
      stop("Item `", col, "` contains values other than 0 and 1. ",
           "iDIFr requires dichotomously scored items.", call. = FALSE)
    }
  }

  cols
}

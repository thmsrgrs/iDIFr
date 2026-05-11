#' Check group structure and cell sizes before running DIF analysis
#'
#' @description
#' Provides a concise summary of the group structure defined by your demographic
#' variables. Reports how many groups meet the recommended minimum cell size,
#' optionally checks which levels of specified variables are fully crossed, and
#' points to `group_details()` and `cross_details()` for full breakdowns.
#'
#' @param data A data frame containing demographic variables.
#' @param group A one-sided formula specifying the grouping variable(s),
#'   using the same syntax as `idifr()`. Example: `~ gender * nationality`.
#' @param min_cell_size Minimum recommended group size. Default is 50.
#' @param cross_by Optional character vector of variable name(s) to check for
#'   complete crossing. For each unique value of the named variable(s), the
#'   function checks whether every intersectional cell containing that value
#'   meets `min_cell_size`. Example: `cross_by = "nationality"` reports which
#'   nationalities are fully crossed across all other demographic variables.
#'   Multiple variables can be supplied: `cross_by = c("nationality", "gender")`.
#' @param plot Logical. If `TRUE` (default), prints a heatmap of cell sizes.
#'   Only applies when there are at least two grouping variables.
#'
#' @return An object of class `idifr_groups` (invisibly), which can be passed
#'   to `merge_groups()`, `group_details()`, or `cross_details()`.
#'
#' @examples
#' \dontrun{
#' # Basic group check
#' grp <- check_groups(my_data, group = ~ gender * nationality * age_band)
#'
#' # With crossing check
#' grp <- check_groups(my_data,
#'                     group    = ~ gender * nationality * age_band,
#'                     cross_by = "nationality")
#'
#' # See full breakdowns
#' group_details(grp)
#' cross_details(grp, cross_by = "nationality")
#' }
#'
#' @seealso [group_details()], [cross_details()], [merge_groups()], [idifr()]
#' @export
check_groups <- function(data,
                         group,
                         min_cell_size = 50,
                         cross_by      = NULL,
                         plot          = TRUE) {

  if (!inherits(group, "formula")) {
    stop("`group` must be a formula. Example: group = ~ gender * nationality",
         call. = FALSE)
  }

  groups <- .build_groups(data, group, min_cell_size, verbose = FALSE)

  vars   <- groups$vars
  n_vars <- length(vars)
  tbl    <- groups$cell_table

  # --- Header -----------------------------------------------------------------

  cli::cli_h1("Group structure summary")

  if (n_vars == 1) {
    cli::cli_text("Single demographic variable: {.field {vars}}")
  } else {
    cli::cli_text(
      "Intersectional design: {.field {paste(vars, collapse = ' \u00d7 ')}}"
    )
    cli::cli_text("{groups$n_groups} groups from crossing {n_vars} variables")
  }
  cat("\n")

  # --- Part 1: Overall cell size summary (one line) ---------------------------

  n_small   <- sum(tbl$n < min_cell_size)
  n_total   <- nrow(tbl)
  pct_small <- round(100 * n_small / n_total)

  if (n_small == 0) {
    cli::cli_alert_success(
      "All {n_total} groups meet the recommended minimum of {min_cell_size}."
    )
  } else {
    cli::cli_alert_warning(
      "{n_small} of {n_total} groups ({pct_small}%) are below the \\
       recommended minimum of {min_cell_size}."
    )
  }

  cli::cli_text(
    "Run {.code group_details(grp)} to see the full per-group breakdown."
  )
  if (n_vars > 1 && is.null(cross_by)) {
    cli::cli_text(
      "Tip: use {.code cross_by = '{vars[1]}'} in {.code check_groups()} \\
       to check which levels are fully crossed across all other variables."
    )
  }
  cat("\n")

  # --- Part 2: Crossing check (optional) --------------------------------------

  crossing_results <- NULL

  if (!is.null(cross_by)) {

    bad_vars <- setdiff(cross_by, vars)
    if (length(bad_vars) > 0) {
      cli::cli_alert_warning(
        "These cross_by variable(s) are not in the group formula and will \\
         be ignored: {paste(bad_vars, collapse = ', ')}"
      )
      cross_by <- intersect(cross_by, vars)
    }

    if (length(cross_by) > 0) {

      cli::cli_h2("Complete crossing check")

      crossing_results <- lapply(cross_by, function(cv) {
        .check_crossing(tbl, cv, vars, min_cell_size)
      })
      names(crossing_results) <- cross_by

      for (cv in cross_by) {

        res         <- crossing_results[[cv]]
        n_crossed   <- sum(res$fully_crossed)
        n_levels    <- nrow(res)
        pct_crossed <- round(100 * n_crossed / n_levels)

        cli::cli_text(
          "Variable: {.field {cv}}  \u2014  \\
           {n_crossed} of {n_levels} levels ({pct_crossed}%) fully crossed"
        )

        if (n_crossed == 0) {
          cli::cli_alert_warning(
            "No levels of {.field {cv}} are fully crossed."
          )
        } else {
          crossed_levels <- res$level[res$fully_crossed]
          cli::cli_alert_success(
            "Fully crossed: {paste(as.character(crossed_levels), collapse = ', ')}"
          )
        }

        cli::cli_text(
          "Run {.code cross_details(grp, cross_by = '{cv}')} to see the \\
           full breakdown by {cv} level."
        )
        cat("\n")
      }
    }
  }

  # Store crossing results in the groups object for use by cross_details()
  groups$crossing_results <- crossing_results
  groups$cross_by         <- cross_by

  # --- Plot if requested ------------------------------------------------------

  if (plot && n_vars >= 2) {
    p <- .plot_cell_heatmap(groups, min_cell_size)
    print(p)
  }

  invisible(groups)
}


#' Full per-group cell size breakdown
#'
#' @description
#' Prints a detailed table showing the cell size for every intersectional
#' group, flagging those below the recommended minimum. This is the full
#' breakdown that `check_groups()` summarises in a single line.
#'
#' @param grp An `idifr_groups` object from `check_groups()`.
#' @param min_cell_size Minimum recommended group size. Overrides the stored
#'   value if supplied.
#'
#' @return The `idifr_groups` object, invisibly.
#'
#' @examples
#' \dontrun{
#' grp <- check_groups(my_data, group = ~ gender * nationality * age_band)
#' group_details(grp)
#' }
#'
#' @seealso [check_groups()], [cross_details()]
#' @export
group_details <- function(grp, min_cell_size = NULL) {

  if (!inherits(grp, "idifr_groups")) {
    stop("`grp` must be an `idifr_groups` object from `check_groups()`.",
         call. = FALSE)
  }

  if (is.null(min_cell_size)) min_cell_size <- grp$min_cell_size

  vars <- grp$vars
  tbl  <- grp$cell_table

  cli::cli_h1("Full group breakdown")
  cli::cli_text(
    "Grouping: {.field {paste(vars, collapse = ' \u00d7 ')}}"
  )
  cli::cli_text("Minimum recommended cell size: {min_cell_size}")
  cat("\n")

  # Build group labels
  group_label <- do.call(
    paste,
    c(lapply(vars, function(v) tbl[[v]]), sep = "  \u00d7  ")
  )
  max_label <- max(nchar(group_label))

  # Header row
  cat(formatC("Group", width = max_label + 2, flag = "-"), "Status\n")
  cat(strrep("\u2500", max_label + 40), "\n")

  # One row per group
  for (i in seq_len(nrow(tbl))) {
    lbl <- formatC(group_label[i], width = max_label + 2, flag = "-")
    n_i <- tbl$n[i]
    if (n_i < min_cell_size) {
      status <- paste0(
        "\u26a0  n = ", n_i,
        "  (below recommended minimum of ", min_cell_size, ")"
      )
    } else {
      status <- paste0("\u2713  n = ", n_i)
    }
    cat(lbl, status, "\n")
  }

  cat("\n")

  n_small <- sum(tbl$n < min_cell_size)
  n_total <- nrow(tbl)

  if (n_small == 0) {
    cli::cli_alert_success("All {n_total} groups meet the recommended minimum.")
  } else {
    cli::cli_alert_warning(
      "{n_small} of {n_total} groups are below the recommended minimum."
    )
  }
  
  cli::cli_text(
    "Use {.code merge_groups(grp, var = list('new' = c('old1', 'old2')))} \\
     to merge groups before running your analysis."
  )
  cat("\n")
  
  invisible(grp)
}


#' Full crossing breakdown for a demographic variable
#'
#' @description
#' For each unique level of the specified variable, shows whether every
#' intersectional cell containing that level meets the minimum cell size.
#' One row per level, showing how many cells are adequate and the smallest
#' cell size observed.
#'
#' @param grp An `idifr_groups` object from `check_groups()`.
#' @param cross_by Character vector of variable name(s) to check. Must match
#'   variables in the group formula.
#' @param min_cell_size Minimum recommended group size. Overrides the stored
#'   value if supplied.
#'
#' @return The `idifr_groups` object, invisibly.
#'
#' @examples
#' \dontrun{
#' grp <- check_groups(my_data,
#'                     group    = ~ gender * nationality * age_band,
#'                     cross_by = "nationality")
#' cross_details(grp, cross_by = "nationality")
#' }
#'
#' @seealso [check_groups()], [group_details()]
#' @export
cross_details <- function(grp, cross_by, min_cell_size = NULL) {

  if (!inherits(grp, "idifr_groups")) {
    stop("`grp` must be an `idifr_groups` object from `check_groups()`.",
         call. = FALSE)
  }

  if (missing(cross_by) || is.null(cross_by)) {
    stop(
      "You must supply `cross_by`.\n",
      "  Example: cross_details(grp, cross_by = 'nationality')",
      call. = FALSE
    )
  }

  if (is.null(min_cell_size)) min_cell_size <- grp$min_cell_size

  vars <- grp$vars
  tbl  <- grp$cell_table

  bad_vars <- setdiff(cross_by, vars)
  if (length(bad_vars) > 0) {
    stop(
      "These variables are not in the group formula: ",
      paste(bad_vars, collapse = ", "),
      call. = FALSE
    )
  }

  for (cv in cross_by) {

    res <- .check_crossing(tbl, cv, vars, min_cell_size)

    cli::cli_h1("Crossing details: {.field {cv}}")
    cli::cli_text(
      "One row per level of {.field {cv}}. \\
       'Cells OK' shows how many of that level's intersectional cells \\
       meet the minimum of {min_cell_size}. \\
       'Min n' is the smallest cell for that level."
    )
    cat("\n")

    # Column widths
    max_level <- max(nchar(as.character(res$level)))
    max_level <- max(max_level, nchar("Level"))

    # Header
    cat(
      formatC("Level",    width = max_level + 2, flag = "-"),
      formatC("Status",   width = 22,            flag = "-"),
      formatC("Cells OK", width = 12,            flag = "-"),
      formatC("Min n",    width = 8,             flag = "-"),
      "\n"
    )
    cat(strrep("\u2500", max_level + 48), "\n")

    # One row per level
    for (i in seq_len(nrow(res))) {
      lbl <- formatC(as.character(res$level[i]),
                     width = max_level + 2, flag = "-")
      if (res$fully_crossed[i]) {
        status <- formatC("\u2713  Fully crossed",     width = 22, flag = "-")
      } else {
        status <- formatC("\u26a0  Not fully crossed", width = 22, flag = "-")
      }
      cells_ok <- formatC(
        paste0(res$n_ok[i], " / ", res$n_cells[i]),
        width = 12, flag = "-"
      )
      min_n <- formatC(res$min_n[i], width = 8, flag = "-")

      cat(lbl, status, cells_ok, min_n, "\n")
    }

    cat("\n")

    n_crossed <- sum(res$fully_crossed)
    n_levels  <- nrow(res)

    if (n_crossed == n_levels) {
      cli::cli_alert_success(
        "All {n_levels} levels of {.field {cv}} are fully crossed."
      )
    } else {
      not_crossed <- res$level[!res$fully_crossed]
      cli::cli_alert_warning(
        "{n_crossed} of {n_levels} levels of {.field {cv}} are fully crossed."
      )
      cli::cli_text(
        "Not fully crossed: {paste(as.character(not_crossed), collapse = ', ')}"
      )
    }
    
    cli::cli_text(
      "Use {.code merge_groups(grp, var = list('new' = c('old1', 'old2')))} \\
       to merge groups before running your analysis."
    )
    cat("\n")
  }
  
  invisible(grp)
}


# --- Internal: check crossing for one variable --------------------------------

.check_crossing <- function(tbl, cross_var, all_vars, min_cell_size) {

  levels_cv <- unique(tbl[[cross_var]])

  result <- lapply(levels_cv, function(lv) {
    rows    <- tbl[tbl[[cross_var]] == lv, ]
    n_cells <- nrow(rows)
    n_ok    <- sum(rows$n >= min_cell_size)
    min_n   <- min(rows$n)
    data.frame(
      level         = lv,
      n_cells       = n_cells,
      n_ok          = n_ok,
      min_n         = min_n,
      fully_crossed = (n_ok == n_cells),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, result)
}


#' Merge sparse groups
#'
#' @description
#' Combines sparse intersectional cells by collapsing levels of one or more
#' demographic variables. Returns a modified data frame ready to pass back
#' to `idifr()` or `check_groups()`.
#'
#' @param groups An `idifr_groups` object from `check_groups()`, or a data
#'   frame (in which case `grp_formula` must also be supplied).
#' @param grp_formula A formula, required only if `groups` is a raw data frame.
#' @param ... Named arguments specifying merge rules. Each should be named
#'   after a demographic variable, with a named list mapping new level names
#'   to vectors of old level names.
#' @param min_cell_size Minimum cell size to validate against after merging.
#'
#' @return The original data frame with recoded grouping variable(s).
#'
#' @examples
#' \dontrun{
#' grp <- check_groups(my_data, group = ~ nationality * age_band)
#'
#' merged_data <- merge_groups(
#'   grp,
#'   age_band = list("18-30" = c("18-24", "25-30"))
#' )
#' }
#'
#' @export
merge_groups <- function(groups, grp_formula = NULL, ..., min_cell_size = 50) {

  if (inherits(groups, "idifr_groups")) {
    data        <- groups$data
    grp_formula <- groups$formula
  } else if (is.data.frame(groups)) {
    if (is.null(grp_formula)) {
      stop(
        "When supplying a data frame, you must also supply `grp_formula`.\n",
        "  Example: merge_groups(my_data, grp_formula = ~ gender * nationality, ...)",
        call. = FALSE
      )
    }
    data <- groups
  } else {
    stop("`groups` must be an idifr_groups object or a data frame.",
         call. = FALSE)
  }

  merge_rules <- list(...)

  if (length(merge_rules) == 0) {
    cli::cli_alert_info("No merge rules supplied. Returning data unchanged.")
    cli::cli_text("Supply named arguments to specify merges, e.g.:")
    cli::cli_text(
      "  merge_groups(grp, age_band = list('18-30' = c('18-24', '25-30')))"
    )
    return(invisible(data))
  }

  for (var in names(merge_rules)) {
    if (!var %in% names(data)) {
      stop("Variable '", var, "' not found in data.", call. = FALSE)
    }
    rules   <- merge_rules[[var]]
    new_var <- as.character(data[[var]])
    for (new_level in names(rules)) {
      old_levels <- rules[[new_level]]
      new_var[new_var %in% old_levels] <- new_level
    }
    data[[var]] <- factor(new_var)
    cli::cli_alert_success(
      "Merged {length(unlist(rules))} levels of {.field {var}} into \\
       {length(rules)} new level(s)."
    )
  }

  cli::cli_text("\nUpdated group summary after merging:")
  check_groups(data, group = grp_formula, min_cell_size = min_cell_size,
               plot = FALSE)

  invisible(data)
}


# --- Internal: build groups object --------------------------------------------

.build_groups <- function(data, group, min_cell_size = 50, verbose = TRUE) {

  vars <- all.vars(group)

  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "These grouping variables were not found in `data`: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  # Warn for variables that carry no information (only 1 unique level)
  for (v in vars) {
    n_lev <- length(unique(data[[v]]))
    if (n_lev < 2L) {
      cli::cli_alert_warning(
        "Variable {.field {v}} has only one level and will not contribute \\
         to group differentiation. Consider removing it from the group formula."
      )
    }
  }

  if (length(vars) == 1) {
    group_vec <- as.character(data[[vars]])
  } else {
    group_vec <- do.call(
      paste,
      c(lapply(vars, function(v) as.character(data[[v]])), sep = " \u00d7 ")
    )
  }

  group_factor <- factor(group_vec)
  group_levels <- levels(group_factor)
  n_groups     <- length(group_levels)

  if (n_groups < 2L) {
    stop(
      "The group formula resulted in only 1 group. ",
      "At least 2 groups are required for DIF analysis. ",
      "Check that your grouping variable(s) have more than one unique value.",
      call. = FALSE
    )
  }

  cell_counts <- as.data.frame(table(group_factor), stringsAsFactors = FALSE)
  names(cell_counts) <- c(".group_label", "n")

  cell_table <- unique(data.frame(
    data[vars],
    .group_label = group_vec,
    stringsAsFactors = FALSE
  ))
  cell_table <- merge(cell_table, cell_counts, by = ".group_label", all.x = TRUE)
  cell_table <- cell_table[order(cell_table$.group_label), ]
  rownames(cell_table) <- NULL

  has_small_cells <- any(cell_counts$n < min_cell_size)
  small_cells     <- cell_counts[cell_counts$n < min_cell_size, , drop = FALSE]

  structure(
    list(
      formula          = group,
      vars             = vars,
      group_vector     = group_factor,
      group_levels     = group_levels,
      n_groups         = n_groups,
      cell_table       = cell_table,
      cell_counts      = cell_counts,
      has_small_cells  = has_small_cells,
      small_cells      = small_cells,
      min_cell_size    = min_cell_size,
      crossing_results = NULL,
      cross_by         = NULL,
      data             = data
    ),
    class = "idifr_groups"
  )
}


# --- Internal: cell size heatmap ----------------------------------------------

.plot_cell_heatmap <- function(groups, min_cell_size) {

  vars <- groups$vars
  tbl  <- groups$cell_table

  x_var <- vars[1]
  y_var <- vars[2]

  p <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], fill = n)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = n),
      size = 3.5, colour = "white", fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#c0392b", mid = "#f39c12", high = "#27ae60",
      midpoint = min_cell_size, name = "n"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title    = "Group cell sizes",
      subtitle = paste0("Red = below recommended minimum (n < ", min_cell_size, ")"),
      x = x_var, y = y_var
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text  = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (length(vars) >= 3) {
    p <- p + ggplot2::facet_wrap(~ .data[[vars[3]]], labeller = "label_both")
  }

  p
}


# --- S3 methods ---------------------------------------------------------------

#' @export
print.idifr_groups <- function(x, ...) {
  check_groups(x$data, group = x$formula, min_cell_size = x$min_cell_size,
               cross_by = x$cross_by, plot = FALSE)
  invisible(x)
}

#' @export
plot.idifr_groups <- function(x, ...) {
  .plot_cell_heatmap(x, x$min_cell_size)
}

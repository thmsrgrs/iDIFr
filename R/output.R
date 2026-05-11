# Combine results from multiple DIF methods
# Internal function.

.combine_results <- function(results_list, item_cols) {

  if (length(results_list) == 0) {
    return(list(results = data.frame(), group_direction = NULL))
  }

  item_dfs    <- list()
  dir_tables  <- list()

  for (m in names(results_list)) {
    res <- results_list[[m]]

    # LR and ID return a list; LRT returns a data frame directly
    if (is.list(res) && !is.data.frame(res)) {
      item_dfs[[m]]   <- res$item_results
      if (!is.null(res$group_direction)) {
        dir_tables[[m]] <- res$group_direction
      }
    } else {
      item_dfs[[m]] <- res
    }
  }

  # Stack item-level results (bind_rows fills missing columns with NA)
  combined <- dplyr::bind_rows(lapply(names(item_dfs), function(m) {
    df       <- item_dfs[[m]]
    df$item  <- factor(df$item, levels = item_cols)
    df
  }))
  combined <- combined[order(combined$item, combined$method), ]
  rownames(combined) <- NULL

  # Stack direction tables if any
  group_direction <- if (length(dir_tables) > 0) {
    do.call(rbind, dir_tables)
  } else {
    NULL
  }

  list(results = combined, group_direction = group_direction)
}


#' Print method for idifr objects
#
#' @param x An `idifr` object.
#' @param ... Ignored.
#' @export
print.idifr <- function(x, ...) {

  cli::cli_h1("iDIFr -- Intersectional DIF Analysis")

  cat("\n")
  cat("Methods:  ", paste(x$method, collapse = ", "), "\n")
  cat("Items:    ", length(x$items), "\n")
  cat("Groups:   ", x$groups$n_groups,
      paste0("(", paste(x$groups$vars, collapse = " \u00d7 "), ")"), "\n")
  cat("Design:   ", if (!is.null(x$design)) x$design else "unknown", "\n")
  cat("Alpha:    ", x$alpha,
      if (x$p_adjust != "none") paste0("(adjusted: ", x$p_adjust, ")") else "", "\n")
  cat("\n")

  if (x$groups$has_small_cells) {
    cli::cli_alert_warning(
      "{sum(x$groups$cell_counts$n < x$groups$min_cell_size)} group(s) below \\
       recommended minimum (n < {x$groups$min_cell_size})."
    )
    cat("\n")
  }

  # --- Flagged items summary --------------------------------------------------

  flagged <- x$results[!is.na(x$results$flagged) & x$results$flagged, ]

  if (nrow(flagged) == 0) {
    cli::cli_alert_success("No items flagged for DIF across any method.")
    cat("\n")
    return(invisible(x))
  }

  cli::cli_h2(paste0(
    length(unique(flagged$item)), " item(s) flagged for DIF"
  ))
  cat("\n")

  # Print flagged items, leading with effect size
  for (item_i in unique(as.character(flagged$item))) {

    item_rows <- flagged[as.character(flagged$item) == item_i, ]
    cat(cli::col_red(paste0("\u25cf  ", item_i)), "\n")

    for (j in seq_len(nrow(item_rows))) {
      r <- item_rows[j, ]

      if (r$method == "LR") {
        es_str  <- paste0("delta-R2 = ", r$delta_r2, "  [", r$ets_class, "]")
        dif_str <- if (!is.na(r$dif_type)) r$dif_type else ""

      } else if (r$method == "LRT") {
        es_str  <- paste0("std-chi  = ", r$std_chi, "  [", r$es_class, "]")
        dif_str <- if (!is.na(r$dif_type)) r$dif_type else ""

      } else if (r$method == "ID") {
        es_str  <- paste0("delta-R2 = ", r$delta_r2_omnibus,
                          "  [", r$ets_class_omnibus, "]")
        dif_str <- paste0(
          "Source: ", if (!is.na(r$dif_source)) r$dif_source else "NA",
          "  (main=",         round(r$delta_r2_main,         3),
          "  twoway=",        round(r$delta_r2_twoway,        3),
          "  intersection=",  round(r$delta_r2_intersection,  3), ")"
        )

      } else {
        es_str  <- ""
        dif_str <- ""
      }

      p_str <- if (!is.na(r$p_adj)) {
        paste0("p_adj = ", formatC(r$p_adj, digits = 3, format = "f"))
      } else "p = NA"

      cat(sprintf(
        "    %-8s  %s  %s  %s\n",
        r$method, es_str, p_str, dif_str
      ))
    }

    # Direction / discrimination table for this item
    if (!is.null(x$group_direction)) {
      item_dir <- x$group_direction[x$group_direction$item == item_i, ]

      if (nrow(item_dir) > 0) {

        dif_t <- item_dir$dif_type[1]

        if (dif_t == "ID") {

          # --- ID decomposition table ------------------------------------------
          # Pull omnibus delta-R2 from results (same source as the es_str line)
          id_res_row   <- x$results[
            x$results$method == "ID" & as.character(x$results$item) == item_i, ]
          omnibus_r2   <- if (nrow(id_res_row) == 1) id_res_row$delta_r2_omnibus
                          else NA_real_

          header <- paste0(
            "ID decomposition (omnibus delta-R2 = ",
            round(omnibus_r2, 4), "):"
          )
          cat("\n   ", header, "\n")
          cat(sprintf(
            "    %-30s  %-26s  %-10s  %s\n",
            "Test level", "Result", "chi-sq", "delta-R2"
          ))
          cat("    ", strrep("\u2500", 72), "\n", sep = "")

          for (k in seq_len(nrow(item_dir))) {
            d <- item_dir[k, ]
            cat(sprintf(
              "    %-30s  %-26s  %-10s  %s\n",
              d$group,
              if (!is.na(d$direction)) d$direction else "NA",
              if (!is.na(d$value))    sprintf("%.3f",  d$value)    else "NA",
              if (!is.na(d$baseline)) sprintf("%.4f",  d$baseline) else "NA"
            ))
          }
          cat("\n")

        } else {

          # --- LR group direction table ----------------------------------------
          metric    <- item_dir$metric[1]
          is_dichot <- isTRUE(x$design == "dichotomous")
          ref_group <- if (is_dichot) x$groups$group_levels[1] else NULL

          header <- dplyr::case_when(
            dif_t == "Non-uniform" ~
              "Group discrimination vs M1 baseline:",
            grepl("group mean", metric, ignore.case = TRUE) & is_dichot ~
              paste0("Group item difficulty (reference: ", ref_group,
                     ";  mean beta = ", round(item_dir$baseline[1], 3), "):"),
            grepl("group mean", metric, ignore.case = TRUE) ~
              paste0("Group item difficulty vs cross-group mean (mean beta = ",
                     round(item_dir$baseline[1], 3), "):"),
            is_dichot ~
              paste0("P(correct) at mean ability (reference: ", ref_group,
                     ";  mean P = ", round(item_dir$baseline[1], 3), "):"),
            TRUE ~
              paste0("Group advantage at mean ability (cross-group mean P = ",
                     round(item_dir$baseline[1], 3), "):")
          )

          cat("\n   ", header, "\n")
          cat(sprintf(
            "    %-30s  %-26s  %-8s  %s\n",
            "Group", "Direction", "Value", "Deviation"
          ))
          cat("    ", strrep("\u2500", 72), "\n", sep = "")

          for (k in seq_len(nrow(item_dir))) {
            d <- item_dir[k, ]
            val_str <- if (!is.na(d$value))    sprintf("%.3f",  d$value)    else "NA"
            dev_str <- if (!is.na(d$deviation)) sprintf("%+.3f", d$deviation) else "NA"
            dir_str <- if (!is.na(d$direction)) d$direction                  else "NA"
            cat(sprintf(
              "    %-30s  %-26s  %-8s  %s\n",
              d$group, dir_str, val_str, dev_str
            ))
          }
          cat("\n")
        }
      }
    }

    cat("\n")
  }

  # --- Non-flagged items (brief) ----------------------------------------------

  not_flagged <- setdiff(x$items, unique(as.character(flagged$item)))
  if (length(not_flagged) > 0) {
    cli::cli_text(
      "{length(not_flagged)} item(s) showed no DIF: \\
       {paste(not_flagged, collapse = ', ')}"
    )
  }

  cat("\n")
  cli::cli_text("Use {.code summary()} for full results or {.code tidy()} for a flat data frame.")
  cat("\n")

  invisible(x)
}


#' Summary method for idifr objects
#
#' @param object An `idifr` object.
#' @param ... Ignored.
#' @export
summary.idifr <- function(object, ...) {

  cli::cli_h1("iDIFr Results Summary")
  cat("\n")

  res <- object$results

  # Per-method summary
  for (m in object$method) {
    m_res <- res[res$method == m, ]
    n_flagged <- sum(!is.na(m_res$flagged) & m_res$flagged)

    cli::cli_h2(m)
    cat(sprintf(
      "  Items analysed: %d    Items flagged: %d    Proportion: %.1f%%\n\n",
      nrow(m_res), n_flagged, 100 * n_flagged / nrow(m_res)
    ))

    # Effect size distribution
    if (m == "LR" && "ets_class" %in% names(m_res)) {
      tbl <- table(m_res$ets_class)
      cat("  Effect size distribution (ETS):\n")
      for (cls in names(tbl)) {
        cat(sprintf("    %-20s  %d items\n", cls, tbl[cls]))
      }
    } else if ("es_class" %in% names(m_res)) {
      tbl <- table(m_res$es_class)
      cat("  Effect size distribution:\n")
      for (cls in names(tbl)) {
        cat(sprintf("    %-20s  %d items\n", cls, tbl[cls]))
      }
    }
    cat("\n")
  }

  # Concordance across methods if >1 method
  if (length(object$method) > 1) {
    .print_concordance(object)
  }

  invisible(object)
}


.print_concordance <- function(x) {
  cli::cli_h2("Method concordance")
  cat("\n")

  methods <- x$method
  items   <- x$items

  flag_mat <- sapply(methods, function(m) {
    m_res <- x$results[x$results$method == m, ]
    m_res <- m_res[match(items, as.character(m_res$item)), ]
    ifelse(is.na(m_res$flagged), FALSE, m_res$flagged)
  })
  rownames(flag_mat) <- items

  # Items flagged by all methods
  all_methods  <- rowSums(flag_mat) == length(methods)
  some_methods <- rowSums(flag_mat) > 0 & !all_methods
  no_methods   <- rowSums(flag_mat) == 0

  cat(sprintf("  Flagged by all methods:  %d items\n", sum(all_methods)))
  cat(sprintf("  Flagged by some methods: %d items\n", sum(some_methods)))
  cat(sprintf("  Flagged by no methods:   %d items\n", sum(no_methods)))

  if (sum(some_methods) > 0) {
    cat("\n  Items with method disagreement:\n")
    disagree_items <- items[some_methods]
    for (item_i in disagree_items) {
      flags <- flag_mat[item_i, ]
      flagged_by   <- methods[flags]
      unflagged_by <- methods[!flags]
      cat(sprintf(
        "    %s: flagged by [%s], not by [%s]\n",
        item_i,
        paste(flagged_by, collapse = ", "),
        paste(unflagged_by, collapse = ", ")
      ))
    }
  }
  cat("\n")
}


#' Plot method for idifr objects
#
#' @param x An `idifr` object.
#' @param type Plot type: `"items"` (default, one row per item showing effect
#'   sizes across methods), `"concordance"` (method agreement heatmap), or
#'   `"groups"` (cell size heatmap from the group structure).
#' @param ... Ignored.
#' @export
plot.idifr <- function(x, type = "items", ...) {
  switch(type,
    "items"       = .plot_items(x),
    "concordance" = .plot_concordance(x),
    "groups"      = plot(x$groups),
    stop("Unknown plot type. Choose 'items', 'concordance', or 'groups'.",
         call. = FALSE)
  )
}


.plot_items <- function(x) {

  res <- x$results

  # Standardise effect size column across methods
  res$es <- dplyr::case_when(
    res$method == "LR"  ~ res$delta_r2,
    res$method == "ID"  ~ res$delta_r2_omnibus,
    res$method == "LRT" ~ res$std_chi,
    TRUE                ~ NA_real_
  )

  res$es_label <- dplyr::case_when(
    res$method == "LR"  ~ res$ets_class,
    res$method == "ID"  ~ res$ets_class_omnibus,
    res$method == "LRT" ~ res$es_class,
    TRUE                ~ NA_character_
  )

  res$flagged_f <- ifelse(!is.na(res$flagged) & res$flagged, "Flagged", "Not flagged")

  ggplot2::ggplot(
    res,
    ggplot2::aes(
      x    = method,
      y    = reorder(item, dplyr::desc(item)),
      fill = es
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(!is.na(flagged) & flagged, "\u25cf", "")),
      colour = "white", size = 3
    ) +
    ggplot2::scale_fill_gradient2(
      low      = "#f0f4f8",
      mid      = "#f39c12",
      high     = "#c0392b",
      midpoint = 0.035,
      na.value = "#e0e0e0",
      name     = "Effect size"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title    = "DIF effect sizes by item and method",
      subtitle = "\u25cf = flagged (p_adj < alpha and effect size threshold met)",
      x        = "Method",
      y        = "Item"
    ) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9)
    )
}


.plot_concordance <- function(x) {

  methods <- x$method
  items   <- x$items

  flag_mat <- lapply(methods, function(m) {
    m_res <- x$results[x$results$method == m, ]
    m_res <- m_res[match(items, as.character(m_res$item)), ]
    data.frame(
      item    = items,
      method  = m,
      flagged = ifelse(is.na(m_res$flagged), FALSE, m_res$flagged)
    )
  })

  flag_df <- do.call(rbind, flag_mat)

  ggplot2::ggplot(
    flag_df,
    ggplot2::aes(
      x    = method,
      y    = reorder(item, dplyr::desc(item)),
      fill = flagged
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#c0392b", "FALSE" = "#ecf0f1"),
      labels = c("TRUE" = "Flagged", "FALSE" = "Not flagged"),
      name   = ""
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "DIF detection concordance across methods",
      x     = "Method",
      y     = "Item"
    ) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}


#' Extract tidy DIF results
#
#' @description Generic for extracting tidy data frames from iDIFr objects.
#' @param x An object.
#' @param ... Additional arguments passed to methods.
#' @return A data frame.
#' @export
tidy <- function(x, ...) UseMethod("tidy")

#' Return tidy data frame of DIF results
#
#' @description
#' Returns results as a tidy data frame suitable for use with `dplyr`,
#' `ggplot2`, or for export. Use the `table` argument to choose which
#' table to return.
#
#' @param x An `idifr` object.
#' @param table Which table to return. One of:
#'   \describe{
#'     \item{`"results"`}{(default) One row per item per method. Includes
#'       test statistics, p-values, effect sizes, and DIF classification.}
#'     \item{`"direction"`}{One row per group per flagged item. Shows
#'       direction and magnitude of DIF for each group. Only available
#'       when `method` includes `"LR"`.}
#'   }
#' @param ... Ignored.
#
#' @return A data frame.
#
#' @examples
#' \dontrun{
#' # Item-level results
#' tidy(result)
#' tidy(result, table = "results")
#
#' # Group direction table for flagged items
#' tidy(result, table = "direction")
#' }
#
#' @export
tidy.idifr <- function(x, table = "results", ...) {

  table <- match.arg(table, c("results", "direction"))

  if (table == "results") {
    return(x$results)
  }

  if (table == "direction") {
    if (is.null(x$group_direction)) {
      message(
        "No direction table available. Direction tables are produced for ",
        "flagged items when method = 'LR' or 'ID'."
      )
      return(invisible(NULL))
    }
    return(x$group_direction)
  }
}

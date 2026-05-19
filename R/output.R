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

    # LR returns a list; LRT returns a data frame directly
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
        dif_str <- if (!is.na(r$dif_type)) r$dif_type else ""

        if (!is.na(r$dif_type) && r$dif_type == "Non-uniform") {
          # Non-uniform only: show selected non-uniform effect size
          es_str <- paste0(r$nu_es_label, " = ", r$nu_es, "  [", r$nu_es_class, "]")

        } else if (!is.na(r$dif_type) && r$dif_type == "Uniform and Non-uniform") {
          # Both: two-part display — uniform component + non-uniform component
          es_str <- paste0(
            "delta-R2(uniform) = ", r$delta_r2_uniform,
            "  [", .ets_classify(r$delta_r2_uniform), "]",
            "  |  ",
            r$nu_es_label, " = ", r$nu_es, "  [", r$nu_es_class, "]"
          )

        } else {
          # Uniform or None: show uniform delta-R2
          es_str <- paste0("delta-R2(uniform) = ", r$delta_r2, "  [", r$ets_class, "]")
        }

      } else if (r$method == "LRT") {
        dif_str <- if (!is.na(r$dif_type)) r$dif_type else ""

        if (!is.na(r$dif_type) && r$dif_type == "Non-uniform") {
          # Non-uniform only: show MAPPD (or std_chi_nonuniform if available)
          if (!is.na(r$mappd)) {
            es_str <- paste0("MAPPD = ", r$mappd, "  [", r$mappd_class, "]")
          } else {
            es_str <- paste0("std-chi(non-uniform) = ", r$std_chi_nonuniform,
                             "  [", r$es_class_uniform, "]")
          }

        } else if (!is.na(r$dif_type) && r$dif_type == "Uniform and Non-uniform") {
          # Both: uniform component | MAPPD
          u_class <- if (!is.na(r$es_class_uniform)) r$es_class_uniform else r$es_class
          es_str <- paste0(
            "std-chi(uniform) = ", r$std_chi_uniform, "  [", u_class, "]",
            "  |  ",
            "MAPPD = ", r$mappd, "  [", r$mappd_class, "]"
          )

        } else {
          # Uniform or None: show uniform std-chi
          u_class <- if (!is.na(r$es_class_uniform)) r$es_class_uniform else r$es_class
          es_str  <- paste0("std-chi(uniform) = ", r$std_chi_uniform,
                            "  [", u_class, "]")
        }

      } else if (r$method == "MOB") {
        es_str  <- paste0("std-diff = ", r$std_diff, "  [", r$es_class, "]")
        dif_str <- paste0(
          if (!is.na(r$dif_type) && r$dif_type != "None")
            paste0("[", r$dif_type, "]  ") else "",
          "Source: ", if (!is.na(r$dif_source)) r$dif_source else "NA",
          if (!is.na(r$split_variable))
            paste0("  split: ", r$split_variable, " (depth=", r$split_depth, ")")
          else ""
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

        if (dif_t == "MOB") {

          # --- MOB score-residual table(s) -------------------------------------
          # May have one sub-table (Uniform or Non-uniform) or two stacked
          # (Uniform and Non-uniform), distinguished by the `metric` column.
          mob_res_row <- x$results[
            x$results$method == "MOB" & as.character(x$results$item) == item_i, ]
          split_v <- if (nrow(mob_res_row) == 1) mob_res_row$split_variable else NA

          split_suffix <- if (!is.na(split_v))
            paste0(" (primary split: ", split_v, ")") else ""

          for (met in unique(item_dir$metric)) {
            sub_dir <- item_dir[item_dir$metric == met, ]
            is_nu   <- identical(met, "ability region residuals")

            if (is_nu) {
              # --- Non-uniform: low/high ability pattern table -----------------
              cat("\n   MOB ability-region residuals", split_suffix, ":\n", sep = "")
              cat(sprintf(
                "    %-20s  %-14s  %-14s  %s\n",
                "Group", "Low ability", "High ability", "Discrimination"
              ))
              cat("    ", strrep("\u2500", 72), "\n", sep = "")
              for (k in seq_len(nrow(sub_dir))) {
                d        <- sub_dir[k, ]
                low_str  <- if (!is.na(d$value))    sprintf("%+.4f", d$value)    else "NA"
                hi_str   <- if (!is.na(d$baseline)) sprintf("%+.4f", d$baseline) else "NA"
                pat_str  <- if (!is.na(d$direction)) d$direction else "NA"
                cat(sprintf(
                  "    %-20s  %-14s  %-14s  %s\n",
                  d$group, low_str, hi_str, pat_str
                ))
              }
              cat(paste0(
                "    Note: For non-uniform DIF items the item discriminates more",
                " strongly between ability\n",
                "    levels for one group than another. The group labelled",
                " 'More discriminating' experiences\n",
                "    a steeper relationship between ability and item performance.\n"
              ))
              cat("\n")

            } else {
              # --- Uniform: standard group residual table ----------------------
              cat("\n   MOB group score residuals", split_suffix, ":\n", sep = "")
              cat(sprintf(
                "    %-20s  %-15s  %-15s  %s\n",
                "Group", "Mean residual", "Deviation", "Direction"
              ))
              cat("    ", strrep("\u2500", 72), "\n", sep = "")
              for (k in seq_len(nrow(sub_dir))) {
                d <- sub_dir[k, ]
                cat(sprintf(
                  "    %-20s  %-15s  %-15s  %s\n",
                  d$group,
                  if (!is.na(d$value))     sprintf("%.4f",  d$value)     else "NA",
                  if (!is.na(d$deviation)) sprintf("%+.4f", d$deviation) else "NA",
                  if (!is.na(d$direction)) d$direction                   else "NA"
                ))
              }
              cat("\n")
            }
          }

        } else {

          # --- LR group direction table ----------------------------------------
          # May contain one sub-table (Uniform or Non-uniform) or two
          # (Uniform and Non-uniform), stacked with different dif_type values.
          is_dichot  <- isTRUE(x$design == "dichotomous")
          ref_group  <- if (is_dichot) x$groups$group_levels[1] else NULL
          sub_types  <- unique(item_dir$dif_type)

          for (sub_type in sub_types) {
            sub_dir <- item_dir[item_dir$dif_type == sub_type, ]
            metric  <- sub_dir$metric[1]

            header <- dplyr::case_when(
              sub_type == "Non-uniform" ~
                "Group discrimination vs cross-group mean:",
              grepl("^Group item difficulty vs cross-group mean$",
                    metric, ignore.case = TRUE) ~
                "Group item difficulty vs cross-group mean:",
              grepl("group mean", metric, ignore.case = TRUE) & is_dichot ~
                paste0("Group item difficulty (reference: ", ref_group,
                       ";  mean beta = ", round(sub_dir$baseline[1], 3), "):"),
              grepl("group mean", metric, ignore.case = TRUE) ~
                paste0("Group item difficulty vs cross-group mean (mean beta = ",
                       round(sub_dir$baseline[1], 3), "):"),
              is_dichot ~
                paste0("P(correct) at mean ability (reference: ", ref_group,
                       ";  mean P = ", round(sub_dir$baseline[1], 3), "):"),
              TRUE ~
                paste0("Group advantage at mean ability (cross-group mean P = ",
                       round(sub_dir$baseline[1], 3), "):")
            )

            cat("\n   ", header, "\n")
            cat(sprintf(
              "    %-30s  %-26s  %-8s  %s\n",
              "Group", "Direction", "Value", "Deviation"
            ))
            cat("    ", strrep("\u2500", 72), "\n", sep = "")

            for (k in seq_len(nrow(sub_dir))) {
              d <- sub_dir[k, ]
              val_str <- if (!is.na(d$value))     sprintf("%.3f",  d$value)     else "NA"
              dev_str <- if (!is.na(d$deviation))  sprintf("%+.3f", d$deviation) else "NA"
              dir_str <- if (!is.na(d$direction))  d$direction                   else "NA"
              cat(sprintf(
                "    %-30s  %-26s  %-8s  %s\n",
                d$group, dir_str, val_str, dev_str
              ))
            }
            cat("\n")
          }
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

  # --- ICA section -----------------------------------------------------------

  if (!is.null(x$ica)) {
    cat("\n")
    cli::cli_h2("ICA Classification")

    for (m in x$method) {
      ica_m <- x$ica[x$ica$method == m, ]

      if (length(x$method) > 1) cat(sprintf("\n  Method: %s\n", m))

      .print_ica_block(ica_m)
    }

    cat("\n")
    cli::cli_alert_warning(
      "ICA note: multiple analyses run without cross-analysis p-value \\
       correction. Interpret 'pure_intersection' and 'obscured' \\
       classifications with caution in small samples."
    )
    cat("\n")
  }

  cat("\n")
  cli::cli_text("Use {.code summary()} for full results or {.code tidy()} for a flat data frame.")
  cat("\n")

  invisible(x)
}


.print_ica_block <- function(ica_m) {

  cls_specs <- list(
    list(key = "amplified",         label = "Amplified (DIF in single and intersectional)",
         show_drivers = TRUE),
    list(key = "pure_intersection", label = "Pure intersection (intersectional only)",
         show_drivers = FALSE),
    list(key = "obscured",          label = "Obscured (single-variable only)",
         show_drivers = TRUE),
    list(key = "none",              label = "No DIF",
         show_drivers = FALSE)
  )

  for (spec in cls_specs) {
    rows <- ica_m[ica_m$ica_class == spec$key, ]
    n    <- nrow(rows)
    plural <- if (n == 1) "item" else "items"
    cat(sprintf("\n  %-50s %d %s\n", paste0(spec$label, ":"), n, plural))
    if (n > 0 && spec$show_drivers) {
      for (i in seq_len(n)) {
        driver_str <- if (!is.na(rows$marginal_vars[i]))
          paste0(" (driven by: ", rows$marginal_vars[i], ")")
        else ""
        cat(sprintf("    %s%s\n", rows$item[i], driver_str))
      }
    } else if (n > 0 && spec$key != "none") {
      for (i in seq_len(n)) cat(sprintf("    %s\n", rows$item[i]))
    }
  }
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
    m_res     <- res[res$method == m, ]
    n_flagged <- sum(!is.na(m_res$flagged) & m_res$flagged)

    cli::cli_h2(m)
    cat(sprintf(
      "  Items analysed: %d    Items flagged: %d    Proportion: %.1f%%\n\n",
      nrow(m_res), n_flagged, 100 * n_flagged / nrow(m_res)
    ))

    # DIF type breakdown (flagged items only)
    if ("dif_type" %in% names(m_res) && n_flagged > 0) {
      flagged_type <- m_res[!is.na(m_res$flagged) & m_res$flagged, ]
      valid_types  <- flagged_type$dif_type[!is.na(flagged_type$dif_type)]
      if (length(valid_types) > 0) {
        type_tbl   <- table(valid_types)
        type_order <- c("Uniform", "Non-uniform", "Uniform and Non-uniform")
        cat("  DIF type breakdown:\n")
        for (tp in type_order) {
          if (tp %in% names(type_tbl)) {
            cat(sprintf("    %-25s  %d items\n", tp, type_tbl[[tp]]))
          }
        }
        cat("\n")
      }
    }

    # Effect size distribution (flagged items only)
    if (m == "LR" && "ets_class" %in% names(m_res)) {

      flagged_res  <- m_res[!is.na(m_res$flagged) & m_res$flagged, ]
      uniform_rows <- flagged_res[!is.na(flagged_res$dif_type) &
                                    flagged_res$dif_type %in%
                                      c("Uniform", "Uniform and Non-uniform"), ]
      nu_rows      <- flagged_res[!is.na(flagged_res$dif_type) &
                                    flagged_res$dif_type %in%
                                      c("Non-uniform", "Uniform and Non-uniform"), ]

      if (nrow(uniform_rows) > 0) {
        tbl_u <- table(uniform_rows$ets_class)
        cat("  Effect size distribution for flagged items (ETS, delta-R2):\n")
        for (cls in names(tbl_u)) {
          cat(sprintf("    %-20s  %d items\n", cls, tbl_u[cls]))
        }
      }

      if (nrow(nu_rows) > 0) {
        tbl_nu <- table(nu_rows$nu_es_class)
        cat(sprintf("\n  Non-uniform effect size for flagged items (%s):\n",
                    nu_rows$nu_es_label[1]))
        for (cls in names(tbl_nu)) {
          cat(sprintf("    %-20s  %d items\n", cls, tbl_nu[cls]))
        }
      }

    } else if (m == "LRT" && "es_class_uniform" %in% names(m_res)) {

      flagged_res  <- m_res[!is.na(m_res$flagged) & m_res$flagged, ]
      uniform_rows <- flagged_res[!is.na(flagged_res$dif_type) &
                                    flagged_res$dif_type %in%
                                      c("Uniform", "Uniform and Non-uniform"), ]
      nu_rows      <- flagged_res[!is.na(flagged_res$dif_type) &
                                    flagged_res$dif_type %in%
                                      c("Non-uniform", "Uniform and Non-uniform"), ]

      if (nrow(uniform_rows) > 0) {
        tbl_u <- table(uniform_rows$es_class_uniform)
        cat("  Effect size distribution for flagged items (uniform, std-chi):\n")
        for (cls in names(tbl_u)) {
          cat(sprintf("    %-20s  %d items\n", cls, tbl_u[cls]))
        }
      }

      if (nrow(nu_rows) > 0) {
        tbl_nu <- table(nu_rows$mappd_class)
        cat("\n  Non-uniform effect size for flagged items (MAPPD):\n")
        for (cls in names(tbl_nu)) {
          cat(sprintf("    %-20s  %d items\n", cls, tbl_nu[cls]))
        }
      }

    } else if ("es_class" %in% names(m_res)) {

      flagged_res <- m_res[!is.na(m_res$flagged) & m_res$flagged, ]
      if (nrow(flagged_res) > 0) {
        tbl <- table(flagged_res$es_class)
        cat("  Effect size distribution (flagged items only):\n")
        for (cls in names(tbl)) {
          cat(sprintf("    %-20s  %d items\n", cls, tbl[cls]))
        }
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
    res$method == "LRT" ~ res$std_chi,
    res$method == "MOB"  ~ res$std_diff,
    TRUE                ~ NA_real_
  )

  res$es_label <- dplyr::case_when(
    res$method == "LR"  ~ res$ets_class,
    res$method == "LRT" ~ res$es_class,
    res$method == "MOB"  ~ res$es_class,
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


#' Return tidy data frame of DIF results
#
#' @description
#' Returns results as a tidy data frame suitable for use with `dplyr`,
#' `ggplot2`, or for export. Use the `table` argument to choose which
#' table to return.
#
#' Implements the `tidy` generic from the `generics` package so that
#' `tidy()` works correctly regardless of whether `broom` is also loaded.
#
#' @param x An `idifr` object.
#' @param table Which table to return. One of:
#'   \describe{
#'     \item{`"results"`}{(default) One row per item per method. Includes
#'       test statistics, p-values, effect sizes, and DIF classification.}
#'     \item{`"direction"`}{One row per group per flagged item. Shows
#'       direction and magnitude of DIF for each group. Only available
#'       when `method` includes `"LR"`.}
#'     \item{`"ica"`}{ICA classification table (one row per item per method).
#'       Only available when `idifr()` was called with `ica = TRUE`.}
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
#
#' # ICA classification table (requires ica = TRUE)
#' tidy(result, table = "ica")
#' }
#
#' @importFrom generics tidy
#' @export
tidy.idifr <- function(x, table = "results", ...) {

  table <- match.arg(table, c("results", "direction", "ica"))

  if (table == "results") {
    return(x$results)
  }

  if (table == "direction") {
    if (is.null(x$group_direction)) {
      message(
        "No direction table available. Direction tables are produced for ",
        "flagged items when method = 'LR'."
      )
      return(invisible(NULL))
    }
    return(x$group_direction)
  }

  if (table == "ica") {
    if (is.null(x$ica)) {
      message("No ICA table available. Re-run with ica = TRUE.")
      return(invisible(NULL))
    }
    return(x$ica)
  }
}

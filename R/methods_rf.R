# Random Forest / Structural Change DIF detection (internal)
# Based on Strobl, Wickelmaier & Zeileis (2011) model-based recursive
# partitioning for DIF.
#
# For each item, score residuals are computed (observed - predicted from a
# logistic regression on total score). A recursive CUSUM tree is then built
# over the demographic variables: at each node the variable whose ordering
# produces the most significant structural change (strucchange::sctest) is
# chosen as the split, with Bonferroni correction applied within each node.
#
# @importFrom strucchange efp sctest

.run_rf <- function(data, item_cols, groups, alpha, p_adjust, verbose,
                    irt_scores = FALSE) {

  vars      <- groups$vars
  n_vars    <- length(vars)
  n_items   <- length(item_cols)
  item_data <- data[item_cols]
  min_n     <- max(10L, as.integer(groups$min_cell_size))

  results        <- vector("list", n_items)
  direction_list <- list()

  for (i in seq_along(item_cols)) {

    item_name   <- item_cols[i]
    y           <- item_data[[item_name]]

    other_items <- item_data[, -i, drop = FALSE]
    total_score <- rowSums(other_items, na.rm = TRUE)

    complete    <- !is.na(y)
    y_c         <- y[complete]
    ts_c        <- total_score[complete]

    vars_data <- data[complete, vars, drop = FALSE]
    for (v in vars) vars_data[[v]] <- as.factor(vars_data[[v]])

    # --- Step 1: Score residuals (LR-based) ----------------------------------
    mdat <- data.frame(y_c = y_c, ts_c = ts_c, stringsAsFactors = FALSE)
    fit0 <- tryCatch(
      stats::glm(y_c ~ ts_c, data = mdat, family = stats::binomial()),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(
        stats::glm(y_c ~ ts_c, data = mdat, family = stats::binomial())
      )
    )

    if (is.null(fit0)) {
      results[[i]] <- .rf_na_row(item_name)
      next
    }

    scores <- stats::residuals(fit0, type = "response")  # y - p_hat

    # --- Steps 2-3: Recursive CUSUM tree -------------------------------------
    tree <- .rf_build_tree(
      vars_data  = vars_data,
      scores     = scores,
      vars       = vars,
      alpha_node = alpha,
      min_n      = min_n,
      max_depth  = 3L,
      depth      = 0L
    )

    # --- Step 4: Effect size & classification --------------------------------
    if (!is.na(tree$split_var)) {
      es_info  <- .rf_effect_size(scores, vars_data[[tree$split_var]])
      std_diff <- round(es_info$std_diff, 4)
      es_class <- dplyr::case_when(
        std_diff >= 0.5 ~ "Large",
        std_diff >= 0.2 ~ "Moderate",
        TRUE            ~ "Negligible"
      )
      dif_source <- dplyr::case_when(
        tree$depth >= 3L ~ "intersection",
        tree$depth == 2L ~ "twoway",
        tree$depth == 1L ~ "main",
        TRUE             ~ "none"
      )
    } else {
      std_diff   <- 0
      es_class   <- "Negligible"
      dif_source <- "none"
    }

    results[[i]] <- data.frame(
      item           = item_name,
      method         = "RF",
      split_variable = tree$split_var,
      split_depth    = tree$depth,
      dif_source     = dif_source,
      std_diff       = std_diff,
      es_class       = es_class,
      p_split        = tree$p_split,
      p_overall      = tree$p_split,
      p_adj          = NA_real_,
      flagged        = NA,
      stringsAsFactors = FALSE
    )

    direction_list[[item_name]] <- list(
      scores   = scores,
      split_var = tree$split_var,
      var_data  = if (!is.na(tree$split_var)) vars_data[[tree$split_var]] else NULL,
      depth     = tree$depth
    )

    if (verbose) {
      status <- if (!is.na(tree$split_var)) {
        paste0("[", es_class, "] split=", tree$split_var,
               " depth=", tree$depth)
      } else "[Negligible] No DIF"
      cat(sprintf("  %-20s  %s\n", item_name, status))
    }
  }

  # --- Assemble results and flag --------------------------------------------
  result_df       <- do.call(rbind, results)
  result_df$p_adj <- stats::p.adjust(result_df$p_overall, method = p_adjust)
  result_df$flagged <- !is.na(result_df$p_adj) &
                       result_df$p_adj < alpha  &
                       result_df$es_class != "Negligible"

  # --- Direction tables for flagged items -----------------------------------
  flagged_items   <- result_df$item[!is.na(result_df$flagged) & result_df$flagged]
  group_direction <- NULL

  if (length(flagged_items) > 0) {
    dir_tables <- lapply(as.character(flagged_items), function(nm) {
      stored <- direction_list[[nm]]
      if (is.null(stored) || is.na(stored$split_var) ||
          is.null(stored$var_data)) return(NULL)
      .rf_direction_table(nm, stored$scores, stored$var_data, stored$depth)
    })
    dir_tables <- Filter(Negate(is.null), dir_tables)
    if (length(dir_tables) > 0) {
      group_direction <- do.call(rbind, dir_tables)
      rownames(group_direction) <- NULL
    }
  }

  list(item_results = result_df, group_direction = group_direction)
}


# Recursive CUSUM tree -------------------------------------------------------
# Returns list(split_var, p_split, depth).
# At each node, tests every remaining variable with a Bonferroni-corrected
# OLS-CUSUM test. The best variable (lowest p) becomes the split if it clears
# the corrected threshold. Recursion continues into each level's subgroup.

.rf_build_tree <- function(vars_data, scores, vars, alpha_node,
                           min_n, max_depth, depth = 0L) {

  result <- list(split_var = NA_character_, p_split = NA_real_, depth = depth)

  if (nrow(vars_data) < min_n || depth >= max_depth || length(vars) == 0L) {
    return(result)
  }

  n_vars     <- length(vars)
  alpha_bonf <- alpha_node / n_vars

  p_vals <- setNames(rep(1.0, n_vars), vars)
  for (v in vars) {
    p_vals[v] <- tryCatch({
      ov  <- order(vars_data[[v]])
      so  <- scores[ov]
      efp_obj <- strucchange::efp(so ~ 1, type = "OLS-CUSUM")
      strucchange::sctest(efp_obj)$p.value
    }, error = function(e) 1.0)
  }

  best_var <- vars[which.min(p_vals)]
  best_p   <- p_vals[best_var]

  if (best_p > alpha_bonf) return(result)

  result$split_var <- best_var
  result$p_split   <- best_p
  result$depth     <- depth + 1L

  remaining <- setdiff(vars, best_var)
  if (length(remaining) > 0L) {
    levels_v <- sort(unique(as.character(vars_data[[best_var]])))
    for (lv in levels_v) {
      idx <- as.character(vars_data[[best_var]]) == lv
      if (sum(idx) >= min_n) {
        child <- .rf_build_tree(
          vars_data[idx, , drop = FALSE],
          scores[idx],
          remaining,
          alpha_node,
          min_n, max_depth, depth + 1L
        )
        if (child$depth > result$depth) result$depth <- child$depth
      }
    }
  }

  result
}


# Effect size: max inter-group mean difference / pooled SD -------------------

.rf_effect_size <- function(scores, var_data) {
  lv_chr   <- as.character(var_data)
  levels_v <- sort(unique(lv_chr))

  means  <- sapply(levels_v, function(lv) mean(scores[lv_chr == lv]))
  ns     <- sapply(levels_v, function(lv) sum(lv_chr == lv))

  g_vars <- sapply(levels_v, function(lv) {
    s <- scores[lv_chr == lv]
    if (length(s) > 1L) stats::var(s) else 0
  })

  n_total   <- length(scores)
  n_groups  <- length(levels_v)
  pooled_sd <- if (n_total > n_groups) {
    sqrt(max(sum((ns - 1) * g_vars) / (n_total - n_groups), 1e-10))
  } else sqrt(max(stats::var(scores), 1e-10))

  list(
    std_diff = (max(means) - min(means)) / pooled_sd,
    means    = setNames(means, levels_v)
  )
}


# Direction table for a flagged item -----------------------------------------
# Re-uses the existing group_direction column schema.
# direction: negative residual = Advantaged (item easier than expected);
#            positive residual = Disadvantaged (item harder than expected).

.rf_direction_table <- function(item_name, scores, var_data, depth) {
  lv_chr       <- as.character(var_data)
  levels_v     <- sort(unique(lv_chr))
  overall_mean <- mean(scores)

  means <- sapply(levels_v, function(lv) mean(scores[lv_chr == lv]))

  data.frame(
    item      = item_name,
    dif_type  = "RF",
    group     = levels_v,
    metric    = "mean score residual",
    value     = round(means, 4),
    baseline  = round(overall_mean, 4),
    deviation = round(means - overall_mean, 4),
    direction = ifelse(means < 0, "Advantaged", "Disadvantaged"),
    stringsAsFactors = FALSE
  )
}


# NA placeholder row ----------------------------------------------------------

.rf_na_row <- function(item_name) {
  data.frame(
    item           = item_name,
    method         = "RF",
    split_variable = NA_character_,
    split_depth    = NA_integer_,
    dif_source     = NA_character_,
    std_diff       = NA_real_,
    es_class       = NA_character_,
    p_split        = NA_real_,
    p_overall      = NA_real_,
    p_adj          = NA_real_,
    flagged        = NA,
    stringsAsFactors = FALSE
  )
}

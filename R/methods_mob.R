# MOB -- Model-Based recursive partitioning DIF detection (internal)
# Based on Strobl, Wickelmaier & Zeileis (2011) model-based recursive
# partitioning (MOB) for DIF.
#
# For each item, score residuals are computed (observed - predicted from a
# logistic regression on total score). A recursive CUSUM tree is then built
# over the demographic variables: at each node the variable whose ordering
# produces the most significant structural change (strucchange::sctest) is
# chosen as the split, with Bonferroni correction applied within each node.
#
# @importFrom strucchange efp sctest

.run_mob <- function(data, item_cols, groups, alpha, p_adjust, verbose,
                    irt_scores = FALSE) {

  # Effect-size thresholds — defined once, used for BOTH es_class labels and
  # the flagging criterion so the two are always consistent.
  #
  # std_diff = (max_group_mean - min_group_mean) / pooled_SD on score residuals
  # (a Cohen's d generalised to K groups via the inter-group range).
  #
  # With K groups the expected null value of (max - min)/SD grows as
  #   E[range of K std normals] / sqrt(n/K),
  # so a 0.20 threshold appropriate for two groups is far too permissive when
  # K >= 3 or when large group ability differences (impact) leave residual
  # contamination after total-score conditioning.
  #
  # Calibration against the LR method: LR's ΔR² >= 0.035 gate corresponds
  # to approximately std_diff >= 0.40 for K=3 groups (via η² → Cohen's f
  # → range-based d conversion).  Setting ES_NEGLIGIBLE = 0.35 gives MOB
  # slightly more sensitivity than LR while filtering the noise band that
  # caused over-detection (items with std_diff 0.20–0.35 that are typically
  # impact contamination rather than true item-level DIF).
  ES_NEGLIGIBLE <- 0.35   # std_diff below this → Negligible (not flagged)
  ES_MODERATE   <- 0.70   # std_diff below this → Moderate; above → Large

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
      results[[i]] <- .mob_na_row(item_name)
      next
    }

    scores_raw <- stats::residuals(fit0, type = "response")  # y - p_hat
    scores_abs <- abs(scores_raw)                             # for detection

    # --- Steps 2-3: CUSUM tree (raw + abs scores combined) -------------------
    # Raw residuals detect uniform DIF (mean shift positive→negative across
    # groups) and non-uniform DIF when the ability distribution is asymmetric
    # around the crossing point.  Absolute residuals detect non-uniform DIF
    # more reliably because the high-discrimination group has elevated misfit
    # regardless of direction.  Either metric can miss a DIF type in isolation:
    #   • raw scores fail for moderate-difficulty uniform DIF (E[|y-p̂|] ≈ equal
    #     across groups when p̂ ≈ 0.5, so no mean shift in abs scores either);
    #     actually fail more narrowly when signs cancel.
    #   • abs scores fail for moderate-difficulty uniform DIF for the same reason.
    # Solution: for each split-variable candidate, take the minimum p-value
    # across raw and abs CUSUMs.  The effect-size threshold (std_diff ≥ 0.20)
    # guards against false positives.
    tree <- .mob_build_tree(
      vars_data  = vars_data,
      scores     = scores_raw,
      scores_alt = scores_abs,
      vars       = vars,
      alpha_node = alpha,
      min_n      = min_n,
      max_depth  = 3L,
      depth      = 0L
    )

    # --- Step 4: DIF type via p_raw (level shift?) + crossing pattern --------
    if (!is.na(tree$split_var)) {

      # p_raw: does the raw-score CUSUM fire for this split variable?
      # Fires for uniform DIF (consistent level shift); may not fire for
      # pure crossing-ICC non-uniform DIF (residuals cancel in raw sum).
      p_raw <- tryCatch({
        ov <- order(vars_data[[tree$split_var]])
        sc <- scores_raw[ov]
        if (stats::var(sc) < 1e-10) {
          1.0
        } else {
          strucchange::sctest(
            strucchange::efp(sc ~ 1, type = "OLS-CUSUM")
          )$p.value
        }
      }, error = function(e) {
        # Fallback to Rec-CUSUM if OLS-CUSUM fails numerically.
        tryCatch({
          ov_fb  <- order(vars_data[[tree$split_var]])
          efp_fb <- strucchange::efp(scores_raw[ov_fb] ~ 1, type = "Rec-CUSUM")
          strucchange::sctest(efp_fb)$p.value
        }, error = function(e2) 1.0)
      })

      # crossing: does the logistic discrimination differ between groups?
      # Uses a likelihood-ratio test of ts_c:group interaction in a logistic
      # model.  Pure b-shift (uniform DIF) keeps the same slope → no
      # interaction; a-shift (non-uniform DIF) changes the slope → significant.
      # Returns list(detected, effect) where effect = |β_int| × SD(ts_c).
      crossing_res <- .mob_crossing_detected(
        y_c, ts_c, vars_data[[tree$split_var]]
      )
      crossing <- crossing_res$detected

      # Sequential classification:
      #   p_raw indicates whether a uniform (level-shift) component is present.
      #   crossing indicates whether a non-uniform (crossing-ICC) component
      #   is present.  The two signals are independent.
      dif_type <- dplyr::case_when(
        p_raw < alpha  & crossing  ~ "Uniform and Non-uniform",
        p_raw >= alpha & crossing  ~ "Non-uniform",
        p_raw < alpha  & !crossing ~ "Uniform",
        TRUE                       ~ "None"
      )

      # Effect size: take the largest signal across three metrics.
      # (1) Raw residuals: best for uniform DIF (opposite-sign group means).
      # (2) Abs residuals: captures elevated misfit for the DIF group.
      # (3) For crossing items: max of
      #       (a) Quartile effect — group residual difference at the bottom/top
      #           ability quartile, where the crossing signal is not cancelled.
      #       (b) Standardised interaction coefficient: |β_int| × SD(ts_c),
      #           proportional to discrimination difference and robust even
      #           when the crossing point is far from the ability mean.
      std_diff <- round(max(
        .mob_effect_size(scores_raw, vars_data[[tree$split_var]])$std_diff,
        .mob_effect_size(scores_abs, vars_data[[tree$split_var]])$std_diff,
        if (crossing) max(
          .mob_quartile_effect(scores_raw, vars_data[[tree$split_var]], ts_c),
          crossing_res$effect
        ) else 0
      ), 4)
      es_class <- dplyr::case_when(
        std_diff >= ES_MODERATE   ~ "Large",
        std_diff >= ES_NEGLIGIBLE ~ "Moderate",
        TRUE                      ~ "Negligible"
      )
      dif_source <- dplyr::case_when(
        tree$depth >= 3L ~ "intersection",
        tree$depth == 2L ~ "twoway",
        tree$depth == 1L ~ "main",
        TRUE             ~ "none"
      )
    } else {
      p_raw      <- NA_real_
      crossing   <- FALSE
      dif_type   <- "None"
      std_diff   <- 0
      es_class   <- "Negligible"
      dif_source <- "none"
    }

    results[[i]] <- data.frame(
      item           = item_name,
      method         = "MOB",
      split_variable = tree$split_var,
      split_depth    = tree$depth,
      dif_source     = dif_source,
      dif_type       = dif_type,
      std_diff       = std_diff,
      es_class       = es_class,
      p_split        = tree$p_split,
      p_overall      = tree$p_split,
      p_adj          = NA_real_,
      flagged        = NA,
      stringsAsFactors = FALSE
    )

    direction_list[[item_name]] <- list(
      scores_raw = scores_raw,
      ability    = ts_c,
      split_var  = tree$split_var,
      var_data   = if (!is.na(tree$split_var)) vars_data[[tree$split_var]] else NULL,
      depth      = tree$depth,
      dif_type   = dif_type
    )

    if (verbose) {
      status <- if (!is.na(tree$split_var)) {
        paste0("[", es_class, "] [", dif_type, "] split=", tree$split_var,
               " depth=", tree$depth)
      } else "[Negligible] No DIF"
      cat(sprintf("  %-20s  %s\n", item_name, status))

      # For crossing (non-uniform) items add a per-group discrimination summary.
      if (crossing && !is.na(tree$split_var)) {
        lv_chr   <- as.character(vars_data[[tree$split_var]])
        levels_v <- sort(unique(lv_chr))
        q25_v    <- stats::quantile(ts_c, 0.25, names = FALSE)
        q75_v    <- stats::quantile(ts_c, 0.75, names = FALSE)
        m_lo_vec <- vapply(levels_v, function(lv) {
          idx <- lv_chr == lv
          if (any(idx & ts_c <= q25_v)) mean(scores_raw[idx & ts_c <= q25_v]) else NA_real_
        }, numeric(1L))
        m_hi_vec <- vapply(levels_v, function(lv) {
          idx <- lv_chr == lv
          if (any(idx & ts_c >= q75_v)) mean(scores_raw[idx & ts_c >= q75_v]) else NA_real_
        }, numeric(1L))
        grp_pats <- tolower(.mob_contrast_labels(m_lo_vec, m_hi_vec, levels_v))
        cat(sprintf("    %s\n",
                    paste(paste0(levels_v, ": ", grp_pats), collapse = "  ")))
      }
    }
  }

  # --- Assemble results and flag --------------------------------------------
  result_df       <- do.call(rbind, results)
  result_df$p_adj <- stats::p.adjust(result_df$p_overall, method = p_adjust)
  result_df$flagged <- !is.na(result_df$p_adj) &
                       result_df$p_adj < alpha  &
                       result_df$std_diff >= ES_NEGLIGIBLE

  # --- Direction tables for flagged items -----------------------------------
  flagged_items   <- result_df$item[!is.na(result_df$flagged) & result_df$flagged]
  group_direction <- NULL

  if (length(flagged_items) > 0) {
    dir_tables <- lapply(as.character(flagged_items), function(nm) {
      stored <- direction_list[[nm]]
      if (is.null(stored) || is.na(stored$split_var) ||
          is.null(stored$var_data)) return(NULL)
      .mob_direction_table(nm, stored$scores_raw, stored$var_data, stored$depth,
                           stored$dif_type, stored$ability)
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
# At each node, runs Bonferroni-corrected OLS-CUSUM on raw score residuals
# for every remaining variable; the best variable (lowest p) becomes the
# split if it clears the corrected threshold. DIF TYPE is not determined
# here — the caller uses .mob_crossing_detected() on raw residuals instead.

.mob_build_tree <- function(vars_data, scores, vars, alpha_node,
                            min_n, max_depth, depth = 0L,
                            scores_alt = NULL) {

  result <- list(split_var = NA_character_, p_split = NA_real_, depth = depth)

  if (nrow(vars_data) < min_n || depth >= max_depth || length(vars) == 0L) {
    return(result)
  }

  n_vars     <- length(vars)
  alpha_bonf <- alpha_node / n_vars

  # Helper: OLS-CUSUM p-value for a score vector ordered by variable v.
  # Falls back to Rec-CUSUM on numerical failure; returns 1.0 for near-zero
  # variance or any unrecoverable error.
  .cusum_p <- function(sc_ordered, sc_ordered_orig) {
    tryCatch({
      if (stats::var(sc_ordered) < 1e-10) {
        1.0
      } else {
        efp_obj <- strucchange::efp(sc_ordered ~ 1, type = "OLS-CUSUM")
        strucchange::sctest(efp_obj)$p.value
      }
    }, error = function(e) {
      tryCatch({
        efp_fb <- strucchange::efp(sc_ordered_orig ~ 1, type = "Rec-CUSUM")
        strucchange::sctest(efp_fb)$p.value
      }, error = function(e2) 1.0)
    })
  }

  p_vals <- setNames(rep(1.0, n_vars), vars)
  for (v in vars) {
    ov       <- order(vars_data[[v]])
    sc_prim  <- scores[ov]
    p_prim   <- .cusum_p(sc_prim, sc_prim)

    # When an alternative score vector is supplied (e.g. abs residuals),
    # take the minimum p-value so both uniform and non-uniform DIF are
    # detectable regardless of item difficulty.
    # Apply a 2× Bonferroni correction for running two CUSUM tests (raw and
    # absolute residuals) per variable: min(p1, p2) is stochastically
    # smaller than U[0,1] under H0, so without correction the effective
    # alpha per variable is ~2× the nominal alpha.
    p_alt <- if (!is.null(scores_alt)) {
      sc_alt <- scores_alt[ov]
      .cusum_p(sc_alt, sc_alt)
    } else 1.0

    p_vals[v] <- min(min(p_prim, p_alt) * 2, 1.0)
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
        child <- .mob_build_tree(
          vars_data[idx, , drop = FALSE],
          scores[idx],
          remaining,
          alpha_node,
          min_n, max_depth, depth + 1L,
          scores_alt = if (!is.null(scores_alt)) scores_alt[idx] else NULL
        )
        if (child$depth > result$depth) result$depth <- child$depth
      }
    }
  }

  result
}


# Effect size: max inter-group mean difference / pooled SD -------------------

.mob_effect_size <- function(scores, var_data) {
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


# Quartile effect size for crossing-ICC DIF ------------------------------------
# At the bottom/top ability quartiles the crossing signal is not cancelled by
# residuals from the opposite side of the ICC crossing point.  Returns the
# larger of the two quartile std_diffs.  Called only when crossing = TRUE.

.mob_quartile_effect <- function(scores_raw, var_data, ability) {
  q25 <- stats::quantile(ability, 0.25, names = FALSE)
  q75 <- stats::quantile(ability, 0.75, names = FALSE)

  es_low  <- tryCatch(
    .mob_effect_size(scores_raw[ability <= q25], var_data[ability <= q25]),
    error = function(e) list(std_diff = 0)
  )
  es_high <- tryCatch(
    .mob_effect_size(scores_raw[ability >= q75], var_data[ability >= q75]),
    error = function(e) list(std_diff = 0)
  )

  max(es_low$std_diff, es_high$std_diff, na.rm = TRUE)
}


# Crossing-ICC detection -------------------------------------------------------
# Likelihood-ratio test of the ts_c:group interaction in a logistic model.
#
#   H0: glm(y ~ ts + group, binomial)   -- same discrimination for all groups
#   H1: glm(y ~ ts * group, binomial)   -- discrimination differs by group
#
# A significant interaction means at least one group has a different
# discrimination parameter — the defining feature of crossing-ICC (non-uniform)
# DIF.  Pure b-shift DIF (uniform) keeps the same logit slope and cannot
# trigger this test, unlike linear-residual approximations which also fire for
# any curved residual pattern caused by a b-shift.
#
# alpha_crossing = 0.10 (lenient: second-stage test on already-detected items).

.mob_crossing_detected <- function(y_c, ts_c, var_data, alpha_crossing = 0.10) {
  grp_f <- as.factor(var_data)
  mdat  <- data.frame(y = y_c, ts = ts_c, group = grp_f,
                      stringsAsFactors = FALSE)

  fit_main <- tryCatch(
    stats::glm(y ~ ts + group, data = mdat, family = stats::binomial()),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      stats::glm(y ~ ts + group, data = mdat, family = stats::binomial())
    )
  )
  fit_int  <- tryCatch(
    stats::glm(y ~ ts * group, data = mdat, family = stats::binomial()),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      stats::glm(y ~ ts * group, data = mdat, family = stats::binomial())
    )
  )

  no_cross <- list(detected = FALSE, effect = 0)

  if (is.null(fit_main) || is.null(fit_int)) return(no_cross)

  p_int <- tryCatch({
    av <- stats::anova(fit_main, fit_int, test = "Chisq")
    av[["Pr(>Chi)"]][2]
  }, error = function(e) 1.0)

  if (!isTRUE(!is.na(p_int) && p_int < alpha_crossing)) return(no_cross)

  # Standardised interaction effect: |β_int| × SD(ts_c).
  # Proportional to discrimination difference; robust to off-centre crossings
  # because the slope difference is detectable even when b is in the tails.
  int_coefs <- tryCatch({
    cf <- stats::coef(fit_int)
    cf[grepl(":", names(cf), fixed = TRUE)]
  }, error = function(e) numeric(0))

  effect <- if (length(int_coefs) > 0L) {
    max(abs(int_coefs), na.rm = TRUE) * stats::sd(ts_c)
  } else 0

  list(detected = TRUE, effect = effect)
}


# Comparative discrimination labels for ability-region residuals ---------------
# Labels are assigned by comparing groups to each other, so no two groups
# can receive the same label.
#
# contrast_g = mean_high_g - mean_low_g
# The group with the higher contrast (residuals increase more from low to high
# ability) experiences a more discriminating item — ability predicts performance
# more strongly for that group.
#
# 2 groups:  "More discriminating"  /  "Less discriminating"
# 3+ groups: "More discriminating"  /  "Moderately discriminating"  /
#            "Less discriminating"  (ordinal for middle groups)

.mob_contrast_labels <- function(m_lows, m_his, group_names) {
  contrasts <- m_his - m_lows
  n_grps    <- length(group_names)
  labels    <- setNames(rep(NA_character_, n_grps), group_names)

  valid <- which(!is.na(contrasts))
  if (length(valid) == 0L) return(labels)

  n_valid <- length(valid)
  c_valid <- contrasts[valid]
  # rank 1 = highest contrast = more discriminating
  ranks   <- rank(-c_valid, ties.method = "first")

  labels[valid] <- vapply(seq_len(n_valid), function(i) {
    if (n_valid == 2L) {
      if (ranks[i] == 1L) "More discriminating" else "Less discriminating"
    } else {
      if (ranks[i] == 1L)          "More discriminating"
      else if (ranks[i] == n_valid) "Less discriminating"
      else                          "Moderately discriminating"
    }
  }, character(1L))

  labels
}


# Direction table for a flagged item -----------------------------------------
#
# Uniform DIF:
#   One row per group — mean raw residual vs overall mean.
#   metric = "mean score residual"
#   Columns: group | value (mean) | baseline (overall) | deviation | direction
#
# Non-uniform DIF:
#   One row per group — low/high ability means and crossing pattern label.
#   metric = "ability region residuals"
#   Columns: group | value (low-ab mean) | baseline (high-ab mean) |
#            deviation (overall mean, for reference) | direction (pattern label)
#
# Uniform and Non-uniform: both sub-tables stacked.

.mob_direction_table <- function(item_name, scores_raw, var_data, depth,
                                 dif_type = "Uniform", ability = NULL) {
  lv_chr       <- as.character(var_data)
  levels_v     <- sort(unique(lv_chr))
  overall_mean <- mean(scores_raw)

  # --- Uniform component: one row per group -----------------------------------
  uniform_tbl <- NULL
  if (dif_type %in% c("Uniform", "Uniform and Non-uniform")) {
    means <- sapply(levels_v, function(lv) mean(scores_raw[lv_chr == lv]))
    uniform_tbl <- data.frame(
      item      = item_name,
      dif_type  = "MOB",
      group     = levels_v,
      metric    = "mean score residual",
      value     = round(means, 4),
      baseline  = round(overall_mean, 4),
      deviation = round(means - overall_mean, 4),
      direction = dplyr::case_when(
        is.na(means) ~ NA_character_,
        means > 0   ~ "Advantaged",
        means < 0   ~ "Disadvantaged",
        TRUE        ~ "No difference"
      ),
      stringsAsFactors = FALSE
    )
  }

  # --- Non-uniform component: one row per group (low-ab, high-ab, pattern) ---
  # Uses bottom/top quartiles of ability so the display works even when the
  # ICC crossing point is far from the median.
  # Labels are computed comparatively across all groups (two-pass) so no two
  # groups can receive the same label.
  nu_tbl <- NULL
  if (dif_type %in% c("Non-uniform", "Uniform and Non-uniform") &&
      !is.null(ability)) {
    q25 <- stats::quantile(ability, 0.25)
    q75 <- stats::quantile(ability, 0.75)

    # Pass 1: collect per-group low/high ability means
    m_lows_v <- vapply(levels_v, function(lv) {
      sel <- lv_chr == lv & ability <= q25
      if (sum(sel) > 0L) mean(scores_raw[sel]) else NA_real_
    }, numeric(1L))
    m_his_v <- vapply(levels_v, function(lv) {
      sel <- lv_chr == lv & ability >= q75
      if (sum(sel) > 0L) mean(scores_raw[sel]) else NA_real_
    }, numeric(1L))

    # Pass 2: comparative labels — groups ranked by contrast (m_hi - m_low)
    patterns <- .mob_contrast_labels(m_lows_v, m_his_v, levels_v)

    # value    = low-ability mean residual
    # baseline = high-ability mean residual
    # deviation = overall mean (reference)
    # direction = comparative pattern label
    rows <- lapply(levels_v, function(lv) {
      data.frame(
        item      = item_name,
        dif_type  = "MOB",
        group     = lv,
        metric    = "ability region residuals",
        value     = round(m_lows_v[[lv]], 4),
        baseline  = round(m_his_v[[lv]],  4),
        deviation = round(overall_mean, 4),
        direction = patterns[[lv]],
        stringsAsFactors = FALSE
      )
    })
    nu_tbl <- do.call(rbind, rows)
  }

  rbind(uniform_tbl, nu_tbl)
}


# NA placeholder row ----------------------------------------------------------

.mob_na_row <- function(item_name) {
  data.frame(
    item           = item_name,
    method         = "MOB",
    split_variable = NA_character_,
    split_depth    = NA_integer_,
    dif_source     = NA_character_,
    dif_type       = NA_character_,
    std_diff       = NA_real_,
    es_class       = NA_character_,
    p_split        = NA_real_,
    p_overall      = NA_real_,
    p_adj          = NA_real_,
    flagged        = NA,
    stringsAsFactors = FALSE
  )
}

# Logistic Regression DIF detection (internal)
#
# Implements the standard three-step logistic regression procedure for DIF
# detection. Group membership is treated as a categorical predictor, making
# this method suitable for both two-group and multi-group designs.
#
# For each item, three nested models are compared:
# - M1: item ~ total_score                          (no group effect)
# - M2: item ~ total_score + group                  (uniform DIF)
# - M3: item ~ total_score * group                  (non-uniform DIF)
#
# For flagged items, two additional tables are computed:
#
# Uniform DIF — direction table:
#   Predicted probability of correct response at mean ability for each group,
#   compared to the cross-group mean predicted probability. Groups above the
#   mean are advantaged; groups below are disadvantaged.
#
# Non-uniform DIF — discrimination table:
#   Effective discrimination for each group extracted from M3 (total_score
#   coefficient + group-specific interaction), compared to the M1 baseline
#   discrimination (discrimination under no-DIF null hypothesis). Groups
#   with discrimination above baseline are more discriminating than expected;
#   groups below are less discriminating.
#
# Effect sizes use Nagelkerke delta R-squared (three components):
#   delta_r2_uniform     (M1 vs M2) — pure group-difficulty effect
#   delta_r2_interaction (M2 vs M3) — non-uniform component only
#   delta_r2_omnibus     (M1 vs M3) — total group effect
# The primary reported delta_r2 is component-specific (uniform items use
# delta_r2_uniform; mixed items use delta_r2_omnibus).
#
# For non-uniform DIF the selected effect size (nonuniform_es argument) is
# stored in nu_es / nu_es_class / nu_es_label.  MAPPD (Maximum Absolute
# Predicted Probability Difference across groups over a theta grid) is the
# default: it is on the probability scale and unaffected by the floor effects
# that make delta_r2_interaction very small for crossing ICCs.
#
# ETS classification (A/B/C) applies to Nagelkerke delta R-squared:
# - A (negligible): delta R2 < 0.035
# - B (moderate):   0.035 <= delta R2 < 0.070
# - C (large):      delta R2 >= 0.070
# MAPPD classification: Negligible < 0.05 <= Moderate < 0.10 <= Large.

.run_lr <- function(data, item_cols, groups, alpha, p_adjust, verbose,
                    nonuniform_es = "MAPPD") {

  n_items      <- length(item_cols)
  group_vector <- groups$group_vector
  item_data    <- data[item_cols]

  results        <- vector("list", n_items)
  direction_list <- list()  # one entry per flagged item

  for (i in seq_along(item_cols)) {

    item_name <- item_cols[i]
    y         <- item_data[[item_name]]

    # Purified total score: exclude the focal item
    other_items <- item_data[, -i, drop = FALSE]
    total_score <- rowSums(other_items, na.rm = TRUE)

    # Remove rows with NA on focal item
    complete  <- !is.na(y)
    y_c       <- y[complete]
    ts_c      <- total_score[complete]
    grp_c     <- droplevels(group_vector[complete])

    # Standardise total score for numerical stability
    ts_scaled <- scale(ts_c)[, 1]

    # Fit three models
    m1 <- tryCatch(
      stats::glm(y_c ~ ts_scaled,           family = stats::binomial()),
      error = function(e) NULL
    )
    m2 <- tryCatch(
      stats::glm(y_c ~ ts_scaled + grp_c,   family = stats::binomial()),
      error = function(e) NULL
    )
    m3 <- tryCatch(
      stats::glm(y_c ~ ts_scaled * grp_c,   family = stats::binomial()),
      error = function(e) NULL
    )

    if (is.null(m1) || is.null(m2) || is.null(m3)) {
      results[[i]] <- .lr_na_row(item_name)
      next
    }

    # Likelihood ratio tests
    lrt_m1_m2 <- tryCatch(stats::anova(m1, m2, test = "LRT"), error = function(e) NULL)
    lrt_m1_m3 <- tryCatch(stats::anova(m1, m3, test = "LRT"), error = function(e) NULL)
    lrt_m2_m3 <- tryCatch(stats::anova(m2, m3, test = "LRT"), error = function(e) NULL)

    if (is.null(lrt_m1_m2) || is.null(lrt_m1_m3) || is.null(lrt_m2_m3)) {
      results[[i]] <- .lr_na_row(item_name)
      next
    }

    p_uniform    <- lrt_m1_m2[2, "Pr(>Chi)"]
    p_nonuniform <- lrt_m1_m3[2, "Pr(>Chi)"]
    p_interaction <- lrt_m2_m3[2, "Pr(>Chi)"]

    # Three Nagelkerke delta R-squared values — one per model step:
    #   uniform:     M1 vs M2 — pure group-difficulty effect
    #   interaction: M2 vs M3 — non-uniform component only (after removing uniform)
    #   omnibus:     M1 vs M3 — total group effect (uniform + non-uniform)
    delta_r2_uniform     <- .nagelkerke_delta(m1, m2)
    delta_r2_interaction <- .nagelkerke_delta(m2, m3)
    delta_r2_omnibus     <- .nagelkerke_delta(m1, m3)

    # Interaction chi-sq = M1-M3 deviance minus M1-M2 deviance.
    chi_sq_interaction <- lrt_m1_m3[2, "Deviance"] - lrt_m1_m2[2, "Deviance"]

    # Four-way classification using M1-M2, M1-M3, and M2-M3 LRTs:
    #   M1 vs M3 = omnibus test (any DIF at all)
    #   M2 vs M3 = interaction-only test (non-uniform beyond uniform)
    #   M1 vs M2 = uniform-only test (group main effect)
    dif_type <- dplyr::case_when(
      is.na(p_uniform) || is.na(p_nonuniform) || is.na(p_interaction) ~ NA_character_,
      p_nonuniform >= alpha                                            ~ "None",
      p_interaction <  alpha & p_uniform <  alpha                     ~ "Uniform and Non-uniform",
      p_interaction <  alpha                                          ~ "Non-uniform",
      TRUE                                                             ~ "Uniform"
    )

    # MAPPD: always computed from M3 for reporting, regardless of nonuniform_es.
    mappd <- tryCatch(.mappd(m3, grp_c), error = function(e) NA_real_)

    # Non-uniform effect size: selected metric (MAPPD, delta_r2, or chi_sq).
    if (nonuniform_es == "MAPPD") {
      nu_es       <- mappd
      nu_es_class <- .mappd_classify(mappd)
      nu_es_label <- "MAPPD"
    } else if (nonuniform_es == "delta_r2") {
      nu_es       <- delta_r2_interaction
      nu_es_class <- .ets_classify(delta_r2_interaction)
      nu_es_label <- "delta-R2(interaction)"
    } else {   # "chi_sq"
      nu_es       <- chi_sq_interaction
      nu_es_class <- dplyr::case_when(
        is.na(chi_sq_interaction)  ~ NA_character_,
        chi_sq_interaction >= 10   ~ "Large",
        chi_sq_interaction >= 5    ~ "Moderate",
        TRUE                       ~ "Negligible"
      )
      nu_es_label <- "chi2(interaction)"
    }

    # Primary reported delta_r2 is type-specific (Nagelkerke scale throughout).
    delta_r2 <- dplyr::case_when(
      is.na(dif_type) | dif_type == "None" ~ delta_r2_omnibus,
      dif_type == "Uniform"                ~ delta_r2_uniform,
      dif_type == "Non-uniform"            ~ delta_r2_interaction,
      TRUE                                 ~ delta_r2_omnibus   # Uniform and Non-uniform
    )

    # p_overall for BH adjustment: use the test most relevant to the DIF type.
    # Non-uniform: M1 vs M3 (omnibus) captures the interaction via M3.
    # All others:  M1 vs M2 (uniform) is the primary test.
    p_overall <- dplyr::case_when(
      is.na(dif_type)             ~ p_nonuniform,
      dif_type == "Non-uniform"   ~ p_nonuniform,
      TRUE                        ~ p_uniform
    )

    # ETS classification for uniform items; nu_es_class for non-uniform.
    ets_class <- dplyr::case_when(
      is.na(dif_type) | dif_type == "None"  ~ "A (negligible)",
      dif_type == "Uniform"                 ~ .ets_classify(delta_r2_uniform),
      dif_type == "Non-uniform"             ~ nu_es_class,
      TRUE                                  ~ .ets_classify(delta_r2_omnibus)
    )

    results[[i]] <- data.frame(
      item                 = item_name,
      method               = "LR",
      chi_sq_uniform       = lrt_m1_m2[2, "Deviance"],
      df_uniform           = lrt_m1_m2[2, "Df"],
      p_uniform            = p_uniform,
      chi_sq_nonuniform    = lrt_m1_m3[2, "Deviance"],
      df_nonuniform        = lrt_m1_m3[2, "Df"],
      p_nonuniform         = p_nonuniform,
      chi_sq_interaction   = chi_sq_interaction,
      p_interaction        = p_interaction,
      p_overall            = p_overall,
      delta_r2_uniform     = round(delta_r2_uniform,     4),
      delta_r2_interaction = round(delta_r2_interaction, 4),
      delta_r2_omnibus     = round(delta_r2_omnibus,     4),
      delta_r2             = round(delta_r2,             4),
      mappd                = round(mappd,                4),
      nu_es                = round(nu_es,                4),
      nu_es_class          = nu_es_class,
      nu_es_label          = nu_es_label,
      ets_class            = ets_class,
      dif_type             = dif_type,
      stringsAsFactors     = FALSE
    )

    # --- Direction / discrimination table for this item ----------------------
    # Computed now while models are in scope; stored and attached after
    # p-value adjustment determines which items are truly flagged

    direction_list[[item_name]] <- list(
      dif_type  = dif_type,
      m1        = m1,
      m2        = m2,
      m3        = m3,
      ts_scaled = ts_scaled,
      grp_c     = grp_c,
      alpha     = alpha
    )

    if (verbose) {
      nu_thresh  <- switch(nonuniform_es, MAPPD = 0.05, delta_r2 = 0.035, chi_sq = 3.84)
      has_effect <- delta_r2_uniform >= 0.035 ||
        (!is.na(nu_es) && nu_es >= nu_thresh)
      status <- if (!is.na(p_overall) && p_overall < alpha && has_effect) {
        paste0("[", ets_class, "] ", dif_type, " DIF")
      } else {
        "[A] No DIF"
      }
      cat(sprintf("  %-20s  %s\n", item_name, status))
    }
  }

  # --- P-value adjustment and flagging ---------------------------------------

  result_df       <- do.call(rbind, results)
  result_df$p_adj <- stats::p.adjust(result_df$p_overall, method = p_adjust)

  # Three-part flagging:
  #
  # (a) Uniform — BH-adjusted p + delta_r2_uniform >= 0.035.
  #
  # (b) Uniform-and-Non-uniform — BH-adjusted p + delta_r2_omnibus >= 0.035.
  #     Components are individually smaller; omnibus avoids under-flagging.
  #
  # (c) Non-uniform supplement — bypasses BH because the M2→M3 Nagelkerke
  #     increment and related metrics are inherently small for crossing ICC DIF
  #     even with strong effects. Compensates with a stricter omnibus gate
  #     (p_nonuniform < alpha/2). Threshold is metric-specific.

  nu_threshold <- switch(nonuniform_es, MAPPD = 0.05, delta_r2 = 0.035, chi_sq = 3.84)

  primary <- result_df$p_adj < alpha & (
    (!is.na(result_df$delta_r2_uniform) &
       result_df$delta_r2_uniform >= 0.035) |
    (!is.na(result_df$delta_r2_omnibus) &
       result_df$delta_r2_omnibus >= 0.035 &
       !is.na(result_df$dif_type) &
       result_df$dif_type == "Uniform and Non-uniform")
  )

  nu_supplement <- !is.na(result_df$p_nonuniform) &
    result_df$p_nonuniform < alpha / 2 &
    !is.na(result_df$p_interaction) &
    result_df$p_interaction < alpha &
    !is.na(result_df$nu_es) &
    result_df$nu_es >= nu_threshold

  result_df$flagged <- primary | nu_supplement

  # --- Build direction tables for flagged items only -------------------------

  flagged_items <- result_df$item[!is.na(result_df$flagged) & result_df$flagged]

  group_direction <- NULL

  if (length(flagged_items) > 0) {

    dir_tables <- lapply(as.character(flagged_items), function(item_name) {

      stored   <- direction_list[[item_name]]
      dif_type <- result_df$dif_type[result_df$item == item_name]

      if (is.na(dif_type) || dif_type == "None") return(NULL)

      if (dif_type == "Uniform") {
        .uniform_direction_table(
          item_name = item_name,
          m2        = stored$m2,
          ts_scaled = stored$ts_scaled,
          grp_c     = stored$grp_c
        )
      } else if (dif_type == "Non-uniform") {
        .nonuniform_discrimination_table(
          item_name = item_name,
          m1        = stored$m1,
          m3        = stored$m3,
          grp_c     = stored$grp_c
        )
      } else {
        # Uniform and Non-uniform — return both tables stacked
        rbind(
          .uniform_direction_table(
            item_name = item_name,
            m2        = stored$m2,
            ts_scaled = stored$ts_scaled,
            grp_c     = stored$grp_c
          ),
          .nonuniform_discrimination_table(
            item_name = item_name,
            m1        = stored$m1,
            m3        = stored$m3,
            grp_c     = stored$grp_c
          )
        )
      }
    })

    dir_tables <- Filter(Negate(is.null), dir_tables)

    if (length(dir_tables) > 0) {
      group_direction <- do.call(rbind, dir_tables)
      rownames(group_direction) <- NULL
    }
  }

  # Return both the item-level results and the group direction table
  list(
    item_results    = result_df,
    group_direction = group_direction
  )
}


# --- Uniform DIF: predicted probability direction table ----------------------
#
# For each group, compute predicted P(correct) at mean ability (ts_scaled = 0).
# Compare each group to the cross-group mean predicted probability.
# Deviation above mean = advantaged; below = disadvantaged.

.uniform_direction_table <- function(item_name, m2, ts_scaled, grp_c) {

  group_levels <- levels(grp_c)

  # Build prediction data: one row per group, ts_scaled = 0 (mean ability)
  pred_data <- data.frame(
    ts_scaled = 0,
    grp_c     = factor(group_levels, levels = group_levels)
  )

  pred_probs <- tryCatch(
    stats::predict(m2, newdata = pred_data, type = "response"),
    error = function(e) rep(NA_real_, length(group_levels))
  )

  mean_prob <- mean(pred_probs, na.rm = TRUE)
  deviation <- pred_probs - mean_prob

  direction <- dplyr::case_when(
    is.na(deviation)  ~ NA_character_,
    deviation >  0.05 ~ "Advantaged (easier)",
    deviation < -0.05 ~ "Disadvantaged (harder)",
    deviation >  0    ~ "Slightly advantaged",
    deviation <  0    ~ "Slightly disadvantaged",
    TRUE              ~ "No difference"
  )

  data.frame(
    item        = item_name,
    dif_type    = "Uniform",
    group       = group_levels,
    metric      = "Predicted P(correct) at mean ability",
    value       = round(pred_probs, 3),
    baseline    = round(mean_prob, 3),
    deviation   = round(deviation, 3),
    direction   = direction,
    stringsAsFactors = FALSE
  )
}


# --- Non-uniform DIF: discrimination deviation table -------------------------
#
# Extract effective discrimination for each group from M3:
#   disc_group = coef(ts_scaled) + coef(ts_scaled:group)
# Reference (baseline) discrimination = coef(ts_scaled) from M1, which
# represents the item's discrimination under the no-DIF null hypothesis.
# Deviation from M1 baseline identifies which groups drive non-uniform DIF.

.nonuniform_discrimination_table <- function(item_name, m1, m3, grp_c) {

  group_levels <- levels(grp_c)
  coefs_m3     <- stats::coef(m3)

  # M3 main slope = discrimination for reference group
  ref_slope <- as.numeric(coefs_m3["ts_scaled"])

  # Build effective discrimination for each group from M3
  disc_values <- sapply(group_levels, function(grp) {
    interaction_term <- paste0("ts_scaled:grp_c", grp)
    if (interaction_term %in% names(coefs_m3)) {
      ref_slope + as.numeric(coefs_m3[interaction_term])
    } else {
      ref_slope  # reference group has no interaction term
    }
  })

  # Use cross-group mean as reference so deviations are
  # symmetric: one group must be above, another below.
  mean_disc <- mean(disc_values, na.rm = TRUE)
  deviation <- disc_values - mean_disc

  direction <- dplyr::case_when(
    is.na(deviation)  ~ NA_character_,
    deviation >  0.10 ~ "More discriminating than average",
    deviation < -0.10 ~ "Less discriminating than average",
    deviation >  0    ~ "Slightly more discriminating",
    deviation <  0    ~ "Slightly less discriminating",
    TRUE              ~ "No difference"
  )

  data.frame(
    item        = item_name,
    dif_type    = "Non-uniform",
    group       = group_levels,
    metric      = "Discrimination vs cross-group mean",
    value       = round(disc_values, 3),
    baseline    = round(mean_disc, 3),
    deviation   = round(deviation, 3),
    direction   = direction,
    stringsAsFactors = FALSE
  )
}


# --- Nagelkerke helpers -------------------------------------------------------

.nagelkerke_delta <- function(m_null, m_full) {
  max(0, .nagelkerke_r2(m_full) - .nagelkerke_r2(m_null))
}

.nagelkerke_r2 <- function(model) {
  n      <- length(stats::residuals(model))
  ll_0   <- as.numeric(stats::logLik(
    stats::glm(model$y ~ 1, family = stats::binomial())
  ))
  ll_1   <- as.numeric(stats::logLik(model))
  r2_cox <- 1 - exp(-(2 / n) * (ll_1 - ll_0))
  r2_max <- 1 - exp((2 / n) * ll_0)
  r2_cox / r2_max
}

.ets_classify <- function(delta_r2) {
  dplyr::case_when(
    is.na(delta_r2)  ~ NA_character_,
    delta_r2 < 0.035 ~ "A (negligible)",
    delta_r2 < 0.070 ~ "B (moderate)",
    TRUE             ~ "C (large)"
  )
}

.mappd <- function(m3, grp_c, n_grid = 100) {
  theta_grid   <- seq(-3, 3, length.out = n_grid)
  group_levels <- levels(grp_c)
  pred_list <- lapply(group_levels, function(g) {
    pred_data <- data.frame(
      ts_scaled = theta_grid,
      grp_c     = factor(g, levels = group_levels)
    )
    stats::predict(m3, newdata = pred_data, type = "response")
  })
  max_diff <- 0
  n_groups <- length(group_levels)
  for (i in seq_len(n_groups)) {
    for (j in seq_len(n_groups)) {
      if (i >= j) next
      diff     <- abs(pred_list[[i]] - pred_list[[j]])
      max_diff <- max(max_diff, max(diff))
    }
  }
  max_diff
}

.mappd_classify <- function(mappd) {
  dplyr::case_when(
    is.na(mappd)    ~ NA_character_,
    mappd >= 0.10   ~ "Large",
    mappd >= 0.05   ~ "Moderate",
    TRUE            ~ "Negligible"
  )
}

.lr_na_row <- function(item_name) {
  data.frame(
    item                 = item_name,
    method               = "LR",
    chi_sq_uniform       = NA_real_,
    df_uniform           = NA_integer_,
    p_uniform            = NA_real_,
    chi_sq_nonuniform    = NA_real_,
    df_nonuniform        = NA_integer_,
    p_nonuniform         = NA_real_,
    chi_sq_interaction   = NA_real_,
    p_interaction        = NA_real_,
    p_overall            = NA_real_,
    delta_r2_uniform     = NA_real_,
    delta_r2_interaction = NA_real_,
    delta_r2_omnibus     = NA_real_,
    delta_r2             = NA_real_,
    mappd                = NA_real_,
    nu_es                = NA_real_,
    nu_es_class          = NA_character_,
    nu_es_label          = NA_character_,
    ets_class            = NA_character_,
    dif_type             = NA_character_,
    p_adj                = NA_real_,
    flagged              = NA,
    stringsAsFactors     = FALSE
  )
}

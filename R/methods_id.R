# Intersectional Decomposition (ID) DIF detection (internal)
# Decomposes DIF into: main effects, two-way interactions, and highest-order
# intersection via nested logistic regression LRT comparisons.
#
# Model hierarchy:
#   M0: y ~ total_score                          (baseline)
#   M1: y ~ total_score + main_effects           (adds main effects)
#   M2: y ~ total_score + main + two-way         (adds two-way interactions)
#   M3: y ~ total_score + full_factorial         (adds highest-order intersection)
#
# For single-variable designs, only M0 vs M1 is computed (equivalent to LR).

.run_id <- function(data, item_cols, groups, alpha, p_adjust, verbose) {

  vars          <- groups$vars
  n_vars        <- length(vars)
  n_items       <- length(item_cols)
  item_data     <- data[item_cols]
  is_single_var <- n_vars == 1

  if (is_single_var && verbose) {
    cli::cli_alert_info(
      "ID is equivalent to LR for single-variable designs. \\
       Running main-effect test only."
    )
  }

  results        <- vector("list", n_items)
  direction_list <- list()

  for (i in seq_along(item_cols)) {

    item_name <- item_cols[i]
    y         <- item_data[[item_name]]

    other_items <- item_data[, -i, drop = FALSE]
    total_score <- rowSums(other_items, na.rm = TRUE)

    complete  <- !is.na(y)
    y_c       <- y[complete]
    ts_c      <- total_score[complete]
    ts_scaled <- scale(ts_c)[, 1]

    vars_data <- data[complete, vars, drop = FALSE]
    for (v in vars) vars_data[[v]] <- as.factor(vars_data[[v]])

    mdat <- data.frame(y_c = y_c, ts_scaled = ts_scaled,
                       vars_data, stringsAsFactors = FALSE)

    main_str <- paste(vars, collapse = " + ")
    f0 <- stats::as.formula("y_c ~ ts_scaled")
    f1 <- stats::as.formula(paste("y_c ~ ts_scaled +", main_str))

    if (!is_single_var) {
      pairs      <- combn(vars, 2,
                          FUN = function(x) paste(x, collapse = ":"),
                          simplify = TRUE)
      twoway_str <- paste(pairs, collapse = " + ")
      f2 <- stats::as.formula(
        paste("y_c ~ ts_scaled +", main_str, "+", twoway_str))
      full_str <- paste(vars, collapse = " * ")
      f3 <- stats::as.formula(paste("y_c ~ ts_scaled +", full_str))
    }

    m0 <- tryCatch(
      stats::glm(f0, data = mdat, family = stats::binomial()),
      error = function(e) NULL
    )
    m1 <- tryCatch(
      stats::glm(f1, data = mdat, family = stats::binomial()),
      error = function(e) NULL
    )

    # ---- Single-variable: M0 vs M1 only ---------------------------------------

    if (is_single_var) {

      if (is.null(m0) || is.null(m1)) {
        results[[i]] <- .id_na_row(item_name)
        next
      }

      lrt_main <- tryCatch(
        stats::anova(m0, m1, test = "LRT"),
        error = function(e) NULL
      )
      if (is.null(lrt_main)) {
        results[[i]] <- .id_na_row(item_name)
        next
      }

      p_main   <- lrt_main[2, "Pr(>Chi)"]
      chi_main <- lrt_main[2, "Deviance"]
      df_main  <- lrt_main[2, "Df"]
      dr2_main <- .nagelkerke_delta(m0, m1)
      ets_class  <- .ets_classify(dr2_main)
      dif_source <- if (!is.na(p_main) && p_main < alpha) "main" else "none"

      results[[i]] <- data.frame(
        item                  = item_name,
        method                = "ID",
        chi_sq_main           = round(chi_main, 3),
        df_main               = df_main,
        p_main                = p_main,
        chi_sq_twoway         = NA_real_,
        df_twoway             = NA_integer_,
        p_twoway              = NA_real_,
        chi_sq_intersection   = NA_real_,
        df_intersection       = NA_integer_,
        p_intersection        = NA_real_,
        chi_sq_omnibus        = round(chi_main, 3),
        df_omnibus            = df_main,
        p_omnibus             = p_main,
        delta_r2_main         = round(dr2_main, 4),
        delta_r2_twoway       = NA_real_,
        delta_r2_intersection = NA_real_,
        delta_r2_omnibus      = round(dr2_main, 4),
        ets_class_omnibus     = ets_class,
        dif_source            = dif_source,
        p_overall             = p_main,
        p_adj                 = NA_real_,
        flagged               = NA,
        stringsAsFactors = FALSE
      )

    # ---- Multi-variable: full M0–M3 hierarchy ---------------------------------

    } else {

      m2 <- tryCatch(
        stats::glm(f2, data = mdat, family = stats::binomial()),
        error = function(e) NULL
      )
      m3 <- tryCatch(
        stats::glm(f3, data = mdat, family = stats::binomial()),
        error = function(e) NULL
      )

      if (any(sapply(list(m0, m1, m2, m3), is.null))) {
        results[[i]] <- .id_na_row(item_name)
        next
      }

      lrt_main    <- tryCatch(stats::anova(m0, m1, test = "LRT"),
                              error = function(e) NULL)
      lrt_twoway  <- tryCatch(stats::anova(m1, m2, test = "LRT"),
                              error = function(e) NULL)
      lrt_inter   <- tryCatch(stats::anova(m2, m3, test = "LRT"),
                              error = function(e) NULL)
      lrt_omnibus <- tryCatch(stats::anova(m0, m3, test = "LRT"),
                              error = function(e) NULL)

      if (any(sapply(list(lrt_main, lrt_twoway, lrt_inter, lrt_omnibus),
                     is.null))) {
        results[[i]] <- .id_na_row(item_name)
        next
      }

      p_main    <- lrt_main[2,    "Pr(>Chi)"]
      p_twoway  <- lrt_twoway[2,  "Pr(>Chi)"]
      p_inter   <- lrt_inter[2,   "Pr(>Chi)"]
      p_omnibus <- lrt_omnibus[2, "Pr(>Chi)"]

      chi_main    <- lrt_main[2,    "Deviance"]
      chi_twoway  <- lrt_twoway[2,  "Deviance"]
      chi_inter   <- lrt_inter[2,   "Deviance"]
      chi_omnibus <- lrt_omnibus[2, "Deviance"]

      df_main    <- lrt_main[2,    "Df"]
      df_twoway  <- lrt_twoway[2,  "Df"]
      df_inter   <- lrt_inter[2,   "Df"]
      df_omnibus <- lrt_omnibus[2, "Df"]

      dr2_main    <- .nagelkerke_delta(m0, m1)
      dr2_twoway  <- .nagelkerke_delta(m1, m2)
      dr2_inter   <- .nagelkerke_delta(m2, m3)
      dr2_omnibus <- .nagelkerke_delta(m0, m3)

      ets_class <- .ets_classify(dr2_omnibus)

      dif_source <- dplyr::case_when(
        is.na(p_inter) | is.na(p_twoway) | is.na(p_main) ~ NA_character_,
        p_inter  < alpha ~ "intersection",
        p_twoway < alpha ~ "twoway",
        p_main   < alpha ~ "main",
        TRUE             ~ "none"
      )

      results[[i]] <- data.frame(
        item                  = item_name,
        method                = "ID",
        chi_sq_main           = round(chi_main,   3),
        df_main               = df_main,
        p_main                = p_main,
        chi_sq_twoway         = round(chi_twoway, 3),
        df_twoway             = df_twoway,
        p_twoway              = p_twoway,
        chi_sq_intersection   = round(chi_inter,  3),
        df_intersection       = df_inter,
        p_intersection        = p_inter,
        chi_sq_omnibus        = round(chi_omnibus, 3),
        df_omnibus            = df_omnibus,
        p_omnibus             = p_omnibus,
        delta_r2_main         = round(dr2_main,    4),
        delta_r2_twoway       = round(dr2_twoway,  4),
        delta_r2_intersection = round(dr2_inter,   4),
        delta_r2_omnibus      = round(dr2_omnibus, 4),
        ets_class_omnibus     = ets_class,
        dif_source            = dif_source,
        p_overall             = p_omnibus,
        p_adj                 = NA_real_,
        flagged               = NA,
        stringsAsFactors = FALSE
      )

      # Store decomposition values for direction table (built after flagging)
      direction_list[[item_name]] <- list(
        dr2_main    = dr2_main,
        dr2_twoway  = dr2_twoway,
        dr2_inter   = dr2_inter,
        dr2_omnibus = dr2_omnibus,
        p_main      = p_main,
        p_twoway    = p_twoway,
        p_inter     = p_inter,
        chi_main    = chi_main,
        chi_twoway  = chi_twoway,
        chi_inter   = chi_inter,
        df_main     = df_main,
        df_twoway   = df_twoway,
        df_inter    = df_inter,
        alpha       = alpha
      )
    }

    if (verbose) {
      p_show   <- if (!is_single_var) p_omnibus else p_main
      dr2_show <- if (!is_single_var) dr2_omnibus else dr2_main
      status <- if (!is.na(p_show) && p_show < alpha && dr2_show >= 0.035) {
        paste0("[", .ets_classify(dr2_show), "] dif_source=", dif_source)
      } else {
        "[A] No DIF"
      }
      cat(sprintf("  %-20s  %s\n", item_name, status))
    }
  }

  # --- Assemble result data frame ---------------------------------------------

  result_df       <- do.call(rbind, results)
  result_df$p_adj <- stats::p.adjust(result_df$p_overall, method = p_adjust)
  result_df$flagged <- !is.na(result_df$p_adj) &
                       result_df$p_adj < alpha &
                       !is.na(result_df$delta_r2_omnibus) &
                       result_df$delta_r2_omnibus >= 0.035

  # --- Build group_direction table for flagged multi-variable items -----------

  flagged_items   <- result_df$item[!is.na(result_df$flagged) & result_df$flagged]
  group_direction <- NULL

  if (length(flagged_items) > 0 && !is_single_var) {
    dir_tables <- lapply(as.character(flagged_items), function(nm) {
      stored <- direction_list[[nm]]
      if (is.null(stored)) return(NULL)
      .id_decomposition_table(nm, stored)
    })
    dir_tables <- Filter(Negate(is.null), dir_tables)
    if (length(dir_tables) > 0) {
      group_direction <- do.call(rbind, dir_tables)
      rownames(group_direction) <- NULL
    }
  }

  list(item_results = result_df, group_direction = group_direction)
}


# Build the decomposition summary table for a single flagged item.
# Re-uses the group_direction column schema:
#   group     -> test level label
#   value     -> chi-square statistic
#   baseline  -> delta R2 for that level
#   deviation -> p-value for that level
#   direction -> "Significant" / "Not significant"

.id_decomposition_table <- function(item_name, stored) {
  p_vals  <- c(stored$p_main, stored$p_twoway, stored$p_inter)
  data.frame(
    item      = item_name,
    dif_type  = "ID",
    group     = c("Main effects", "Two-way interactions", "Intersection"),
    metric    = "LRT decomposition",
    value     = round(c(stored$chi_main,   stored$chi_twoway,  stored$chi_inter),  3),
    baseline  = round(c(stored$dr2_main,   stored$dr2_twoway,  stored$dr2_inter),  4),
    deviation = round(p_vals, 4),
    direction = ifelse(
      !is.na(p_vals) & p_vals < stored$alpha,
      "Significant", "Not significant"
    ),
    stringsAsFactors = FALSE
  )
}


# NA placeholder row — returned when model fitting fails for an item.
.id_na_row <- function(item_name) {
  data.frame(
    item                  = item_name,
    method                = "ID",
    chi_sq_main           = NA_real_,
    df_main               = NA_integer_,
    p_main                = NA_real_,
    chi_sq_twoway         = NA_real_,
    df_twoway             = NA_integer_,
    p_twoway              = NA_real_,
    chi_sq_intersection   = NA_real_,
    df_intersection       = NA_integer_,
    p_intersection        = NA_real_,
    chi_sq_omnibus        = NA_real_,
    df_omnibus            = NA_integer_,
    p_omnibus             = NA_real_,
    delta_r2_main         = NA_real_,
    delta_r2_twoway       = NA_real_,
    delta_r2_intersection = NA_real_,
    delta_r2_omnibus      = NA_real_,
    ets_class_omnibus     = NA_character_,
    dif_source            = NA_character_,
    p_overall             = NA_real_,
    p_adj                 = NA_real_,
    flagged               = NA,
    stringsAsFactors = FALSE
  )
}

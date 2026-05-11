# IRT Likelihood Ratio Test DIF detection (internal)
#
# Implements the IRT-LRT approach using the package's own fit_2pl() function --
# no external IRT package required.
#
# Algorithm:
#
# 1. Fit a constrained multigroup 2PL model on the anchor item set
#    (constrain = "items": a and b equal across groups, group ability free).
#    On pass 1 all items are anchors.
#
# 2. For every item j, using the anchor-based posterior as the ability
#    estimate:
#    a. Fit item j's parameters constrained equal across groups (NR, pooled).
#    b. Fit item j's parameters free per group (NR, per group).
#    c. LRT_j = -2 * (LL_j_constrained - LL_j_free), chi-sq with
#       df = 2*(n_groups - 1).
#
# 3. Flag items with p_adj < alpha AND effect size >= moderate.
#
# 4. If purify = TRUE, update the anchor set to exclude flagged items and
#    repeat from step 1 until the flagged set is stable (max max_purify
#    iterations). Using anchor-based posteriors for both constrained and
#    free LLs gives unbiased ability estimates and controls FPR.
#
# Requires fit_2pl() from irt_2pl.R and the Rcpp kernels in src/em_2pl.cpp.

.run_lrt <- function(data, item_cols, groups, anchor, alpha, p_adjust, verbose,
                     purify = TRUE, max_purify = 5) {

  n_items      <- length(item_cols)
  resp_matrix  <- as.matrix(data[item_cols])
  storage.mode(resp_matrix) <- "integer"
  group_vector <- as.character(groups$group_vector)
  group_levels <- groups$group_levels
  n_groups     <- groups$n_groups

  if (n_groups < 2) {
    cli::cli_alert_warning("LRT requires at least 2 groups. Returning NA results.")
    return(.lrt_na_results(item_cols))
  }

  # Validate anchor items
  anchor_cols <- if (!is.null(anchor)) {
    bad <- setdiff(anchor, item_cols)
    if (length(bad) > 0)
      cli::cli_alert_warning("Anchor items not in item list ignored: {paste(bad, collapse=', ')}")
    intersect(anchor, item_cols)
  } else {
    item_cols
  }

  if (length(anchor_cols) < 2) {
    cli::cli_alert_warning("Too few anchor items for LRT. Need at least 2.")
    return(.lrt_na_results(item_cols))
  }

  # ---------------------------------------------------------------------------
  # Iterative purification loop
  # ---------------------------------------------------------------------------

  flagged_prev  <- character(0)
  item_chi_sq   <- rep(NA_real_, n_items)
  item_df       <- 2L * (n_groups - 1L)
  item_p        <- rep(NA_real_, n_items)
  std_chi       <- rep(NA_real_, n_items)
  es_class      <- rep(NA_character_, n_items)
  p_adj_vec     <- rep(NA_real_, n_items)
  flagged_vec   <- rep(NA, n_items)
  group_betas   <- matrix(NA_real_, nrow = n_groups, ncol = n_items)

  for (pass in seq_len(if (purify) max_purify else 1L)) {

    .lrt_log(
      sprintf("Pass %d: constrained model on %d anchor item(s)",
              pass, length(anchor_cols)),
      verbose
    )

    t0 <- proc.time()

    baseline <- tryCatch(
      fit_2pl(
        resp      = resp_matrix[, anchor_cols, drop = FALSE],
        group     = group_vector,
        constrain = "items",
        verbose   = FALSE
      ),
      error = function(e) {
        cli::cli_alert_warning("Constrained model failed: {conditionMessage(e)}")
        NULL
      }
    )

    if (is.null(baseline)) return(.lrt_na_results(item_cols))

    elapsed <- round((proc.time() - t0)["elapsed"], 1)
    .lrt_log(
      sprintf("Constrained model converged (%.1fs)  LL = %.3f",
              elapsed, baseline$loglik),
      verbose, success = TRUE
    )

    # Per-item LLs for ALL items using anchor-based posteriors
    lls <- .lrt_item_lls(
      model        = baseline,
      resp         = resp_matrix,
      group_levels = group_levels,
      n_groups     = n_groups
    )

    # LRT statistics
    item_chi_sq  <- pmax(0, -2 * (lls$constrained - lls$free))
    std_chi      <- sqrt(item_chi_sq / item_df)
    item_p       <- stats::pchisq(item_chi_sq, item_df, lower.tail = FALSE)
    es_class     <- .lrt_es_classify(std_chi, item_df)
    p_adj_vec    <- stats::p.adjust(item_p, method = p_adjust)
    flagged_vec  <- p_adj_vec < alpha & !is.na(es_class) & es_class != "Negligible"
    group_betas  <- lls$group_betas

    flagged_items <- item_cols[!is.na(flagged_vec) & flagged_vec]

    if (!purify || pass == max_purify) break

    # Check convergence
    if (setequal(flagged_items, flagged_prev)) {
      .lrt_log(sprintf("Purification converged after %d pass(es)", pass),
               verbose, success = TRUE)
      break
    }

    flagged_prev <- flagged_items
    new_anchors  <- setdiff(item_cols, flagged_items)

    if (length(new_anchors) < 2) {
      .lrt_log("Too few anchor items remaining -- stopping purification", verbose)
      break
    }

    anchor_cols <- new_anchors
  }

  # ---------------------------------------------------------------------------
  # Build result data frame
  # ---------------------------------------------------------------------------

  result_df <- data.frame(
    item      = item_cols,
    method    = "LRT",
    chi_sq    = round(item_chi_sq, 3),
    df        = item_df,
    p_overall = item_p,
    std_chi   = round(std_chi, 4),
    es_class  = es_class,
    dif_type  = "Uniform",
    stringsAsFactors = FALSE
  )

  result_df$p_adj   <- p_adj_vec
  result_df$flagged <- flagged_vec

  # ---------------------------------------------------------------------------
  # Verbose item-level output
  # ---------------------------------------------------------------------------

  if (verbose) {
    for (i in seq_len(nrow(result_df))) {
      r <- result_df[i, ]
      if (!is.na(r$flagged) && r$flagged) {
        cat(sprintf(
          "  %-20s  \u26a0  chi2=%-8.3f  std_chi=%-6.4f  [%s]\n",
          r$item, r$chi_sq, r$std_chi, r$es_class
        ))
      } else {
        cat(sprintf(
          "  %-20s  \u2713  No DIF  (chi2=%.3f  p=%.3f)\n",
          r$item, r$chi_sq, r$p_adj
        ))
      }
    }
    cat("\n")
  }

  # ---------------------------------------------------------------------------
  # Direction table for flagged items
  # ---------------------------------------------------------------------------

  flagged_names   <- item_cols[!is.na(flagged_vec) & flagged_vec]
  group_direction <- NULL

  if (length(flagged_names) > 0) {

    .lrt_log(
      sprintf("Computing group direction table for %d flagged item(s)...",
              length(flagged_names)),
      verbose
    )

    dir_rows <- lapply(seq_along(flagged_names), function(fi) {
      item_name <- flagged_names[fi]
      j         <- which(item_cols == item_name)
      betas     <- group_betas[, j]
      ref_beta  <- mean(betas, na.rm = TRUE)

      rows <- lapply(seq_along(group_levels), function(gi) {
        g   <- group_levels[gi]
        dev <- betas[gi] - ref_beta
        data.frame(
          item      = item_name,
          dif_type  = "Uniform",
          group     = g,
          metric    = "Item difficulty (beta vs group mean)",
          value     = round(betas[gi], 3),
          baseline  = round(ref_beta, 3),
          deviation = round(dev, 3),
          direction = dplyr::case_when(
            is.na(dev)   ~ NA_character_,
            dev < -0.20  ~ "Advantaged (easier)",
            dev >  0.20  ~ "Disadvantaged (harder)",
            TRUE         ~ "Similar to group mean"
          ),
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, rows)
    })

    group_direction <- do.call(rbind, dir_rows)
    rownames(group_direction) <- NULL
  }

  .lrt_log("LRT analysis complete", verbose, success = TRUE)
  if (verbose) cat("\n")

  list(
    item_results    = result_df,
    group_direction = group_direction
  )
}


# --- Per-item constrained and free log-likelihoods ----------------------------
#
# Uses the anchor-based constrained model's posteriors for both the constrained
# and free LL computation. This gives consistent ability estimates for both
# hypotheses and is the key mechanism behind purification-based FPR control.
#
# For each item j:
#   Constrained LL: a_j, b_j fitted equal across groups (pooled NR)
#   Free LL:        a_jg, b_jg fitted per group (per-group NR)
#
# Both use the SAME anchor-posterior.

.lrt_item_lls <- function(model, resp, group_levels, n_groups) {

  n_items   <- ncol(resp)
  n_nodes   <- length(model$nodes)
  g_int     <- model$g_int           # n_persons, integer group indices
  post      <- model$posterior       # n_persons x n_nodes

  # Group-scaled quadrature nodes (n_groups x n_nodes)
  nodes_g <- outer(sqrt(model$sigma2_vec), model$nodes, "*") +
             matrix(model$mu_vec, nrow = n_groups, ncol = n_nodes)

  # Per-person scaled nodes (n_persons x n_nodes)
  nodes_i <- nodes_g[g_int, , drop = FALSE]

  obs_mat <- !is.na(resp)
  resp_s  <- resp; resp_s[!obs_mat] <- 0L

  ll_c        <- numeric(n_items)
  ll_f        <- numeric(n_items)
  betas_g_mat <- matrix(NA_real_, nrow = n_groups, ncol = n_items)

  for (j in seq_len(n_items)) {

    x_j     <- resp_s[, j]
    obs_j   <- obs_mat[, j]
    obs_idx <- which(obs_j)

    if (length(obs_idx) < 5L) next

    post_obs  <- post[obs_idx, , drop = FALSE]
    x_obs     <- x_j[obs_idx]
    nodes_obs <- nodes_i[obs_idx, , drop = FALSE]
    g_obs     <- g_int[obs_idx]

    # Shared starting values
    pc <- mean(x_obs)
    a0 <- 1.0
    b0 <- qlogis(pmin(pmax(1 - pc, 0.01), 0.99))

    # --- Constrained: a_j, b_j equal across groups ---
    nr_c <- nr_update_pooled(as.double(x_obs), post_obs, nodes_obs, a0, b0)
    a_c  <- .clamp(nr_c[1], 0.1, 5.0)
    b_c  <- .clamp(nr_c[2], -6.0, 6.0)

    eta_c <- a_c * (nodes_obs - b_c)
    P_c   <- pmin(pmax(1 / (1 + exp(-eta_c)), 1e-10), 1 - 1e-10)
    px_c  <- ifelse(matrix(x_obs, length(obs_idx), n_nodes) == 1L, P_c, 1 - P_c)
    ll_c[j] <- sum(log(pmax(rowSums(post_obs * px_c), 1e-300)))

    # --- Free: a_jg, b_jg per group ---
    ll_j_f <- 0

    for (gi in seq_len(n_groups)) {
      g_idx <- which(g_obs == gi)

      if (length(g_idx) < 2L) {
        betas_g_mat[gi, j] <- b_c
        next
      }

      post_gi  <- post_obs[g_idx, , drop = FALSE]
      x_gi     <- x_obs[g_idx]
      nodes_gi <- nodes_obs[g_idx, , drop = FALSE]

      nr_f <- nr_update_pooled(as.double(x_gi), post_gi, nodes_gi, a0, b0)
      a_f  <- .clamp(nr_f[1], 0.1, 5.0)
      b_f  <- .clamp(nr_f[2], -6.0, 6.0)
      betas_g_mat[gi, j] <- b_f

      eta_f <- a_f * (nodes_gi - b_f)
      P_f   <- pmin(pmax(1 / (1 + exp(-eta_f)), 1e-10), 1 - 1e-10)
      px_f  <- ifelse(matrix(x_gi, length(g_idx), n_nodes) == 1L, P_f, 1 - P_f)
      ll_j_f <- ll_j_f + sum(log(pmax(rowSums(post_gi * px_f), 1e-300)))
    }

    ll_f[j] <- ll_j_f
  }

  list(constrained = ll_c, free = ll_f, group_betas = betas_g_mat)
}


# --- Logging helper -----------------------------------------------------------

.lrt_log <- function(msg, verbose, detail = NULL, success = FALSE) {
  if (!verbose) return(invisible(NULL))
  timestamp <- format(Sys.time(), "[%H:%M:%S]")
  symbol    <- if (success) "\u2713" else " "
  cat(sprintf("  %s %s %s", timestamp, symbol, msg))
  if (!is.null(detail)) cat(sprintf("  (%s)", detail))
  cat("\n")
}


# --- Effect size classification -----------------------------------------------
#
# Standardised chi-square index sqrt(X2/df), with df-adjusted thresholds.
# Oshima et al. (1997) thresholds (0.10 moderate, 0.20 large) were derived
# for 2-group designs (df=2). With more groups df increases, so the same
# raw std_chi value has less evidential weight -- thresholds are scaled by
# sqrt(df/2) to maintain equivalent sensitivity across designs.

.lrt_es_classify <- function(std_chi, df = 2) {
  scale_factor    <- sqrt(df / 2)
  threshold_mod   <- 0.10 * scale_factor
  threshold_large <- 0.20 * scale_factor
  dplyr::case_when(
    is.na(std_chi)             ~ NA_character_,
    std_chi >= threshold_large ~ "Large",
    std_chi >= threshold_mod   ~ "Moderate",
    TRUE                       ~ "Negligible"
  )
}


# --- NA helpers ---------------------------------------------------------------

.lrt_na_row <- function(item_name) {
  data.frame(
    item        = item_name,
    method      = "LRT",
    chi_sq      = NA_real_,
    df          = NA_integer_,
    p_overall   = NA_real_,
    std_chi     = NA_real_,
    es_class    = NA_character_,
    dif_type    = NA_character_,
    p_adj       = NA_real_,
    flagged     = NA,
    stringsAsFactors = FALSE
  )
}

.lrt_na_results <- function(item_cols) {
  list(
    item_results    = do.call(rbind, lapply(item_cols, .lrt_na_row)),
    group_direction = NULL
  )
}

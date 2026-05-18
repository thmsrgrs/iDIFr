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
#    estimate, fit three nested models:
#    a. Constrained:       a and b equal across groups  (DIF null)
#    b. Alpha-fixed:       a equal, b free per group    (uniform DIF only)
#    c. Free:              a and b free per group       (omnibus)
#
#    Three LRTs per item:
#      chi_uniform    = -2*(LL_c    - LL_alpha)  df = n_groups - 1
#      chi_nonuniform = -2*(LL_alpha - LL_free)  df = n_groups - 1
#      chi_omnibus    = -2*(LL_c    - LL_free)   df = 2*(n_groups - 1)
#
# 3. BH-correct uniform and non-uniform p-values separately.
#    Classify dif_type; apply type-specific effect-size criterion.
#
# 4. If purify = TRUE, update the anchor set to exclude flagged items and
#    repeat from step 1 until the flagged set is stable (max max_purify
#    iterations).
#
# Requires fit_2pl() from irt_2pl.R and the Rcpp kernels in src/em_2pl.cpp.

.run_lrt <- function(data, item_cols, groups, anchor, alpha, p_adjust, verbose,
                     nonuniform_es = "MAPPD",
                     purify = TRUE, max_purify = 2,
                     cores = NULL) {

  # Resolve nonuniform_es: "delta_r2" is LR-only, fall back to "chi_sq" for LRT
  nonuniform_es <- if (nonuniform_es == "delta_r2") "chi_sq" else nonuniform_es

  # Resolve core count for parallelisation
  n_cores <- max(1L, if (is.null(cores))
    parallel::detectCores(logical = FALSE) - 1L
  else
    as.integer(cores)
  )

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

  flagged_prev <- character(0)
  warm_start   <- NULL  # warm-start params for purification passes

  # Storage for final-pass results (updated each pass)
  lls          <- NULL
  item_df_om   <- 2L * (n_groups - 1L)
  p_adj_om_vec <- rep(NA_real_, n_items)
  flagged_vec  <- rep(NA, n_items)

  for (pass in seq_len(if (purify) max_purify else 1L)) {

    .lrt_log(
      sprintf("Pass %d: constrained model on %d anchor item(s)",
              pass, length(anchor_cols)),
      verbose
    )

    t0 <- proc.time()

    # Build start list aligned to the current anchor item set
    pass_start <- if (!is.null(warm_start) && pass > 1L) {
      keep <- anchor_cols %in% names(warm_start$a)
      if (all(keep)) {
        list(a = warm_start$a[anchor_cols],
             b = warm_start$b[anchor_cols])
      } else NULL
    } else NULL

    baseline <- tryCatch(
      fit_2pl(
        resp      = resp_matrix[, anchor_cols, drop = FALSE],
        group     = group_vector,
        constrain = "items",
        start     = pass_start,
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

    # Capture item parameters for warm-starting the next purification pass.
    # baseline$item_params has columns item, a, b (shared across groups for
    # constrained model).
    warm_start <- list(
      a = stats::setNames(baseline$item_params$a, baseline$item_params$item),
      b = stats::setNames(baseline$item_params$b, baseline$item_params$item)
    )

    # Per-item LLs for ALL items using anchor-based posteriors
    lls <- .lrt_item_lls(
      model        = baseline,
      resp         = resp_matrix,
      group_levels = group_levels,
      n_groups     = n_groups,
      cores        = n_cores
    )

    # Omnibus stats for purification criterion
    chi_om     <- pmax(0, -2 * (lls$constrained - lls$free))
    std_chi_om <- sqrt(chi_om / item_df_om)
    p_om       <- stats::pchisq(chi_om, item_df_om, lower.tail = FALSE)
    p_adj_om_vec <- stats::p.adjust(p_om, method = p_adjust)
    es_class_om  <- .lrt_es_classify(std_chi_om, item_df_om)

    flagged_vec   <- p_adj_om_vec < alpha &
                     !is.na(es_class_om) & es_class_om != "Negligible"
    flagged_items <- item_cols[!is.na(flagged_vec) & flagged_vec]

    if (!purify || pass == max_purify) break

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
  # Three-model decomposition on final lls
  # ---------------------------------------------------------------------------

  df_uniform    <- n_groups - 1L
  df_nonuniform <- n_groups - 1L
  df_omnibus    <- 2L * (n_groups - 1L)

  chi_uniform    <- pmax(0, -2 * (lls$constrained - lls$alpha_fixed))
  chi_nonuniform <- pmax(0, -2 * (lls$alpha_fixed  - lls$free))
  chi_omnibus    <- pmax(0, -2 * (lls$constrained  - lls$free))

  std_chi_uniform    <- sqrt(chi_uniform    / df_uniform)
  std_chi_nonuniform <- sqrt(chi_nonuniform / df_nonuniform)
  std_chi_omnibus    <- sqrt(chi_omnibus    / df_omnibus)

  p_uniform    <- stats::pchisq(chi_uniform,    df_uniform,    lower.tail = FALSE)
  p_nonuniform <- stats::pchisq(chi_nonuniform, df_nonuniform, lower.tail = FALSE)

  p_uniform_adj    <- stats::p.adjust(p_uniform,    method = p_adjust)
  p_nonuniform_adj <- stats::p.adjust(p_nonuniform, method = p_adjust)

  es_class_uniform    <- .lrt_es_classify(std_chi_uniform,    df_uniform)
  es_class_nonuniform <- .lrt_es_classify(std_chi_nonuniform, df_nonuniform)
  es_class_omnibus    <- .lrt_es_classify(std_chi_omnibus,    df_omnibus)

  # MAPPD from free model parameters
  mappd       <- tryCatch(
    .mappd_irt(lls$group_alphas, lls$group_betas, group_levels),
    error = function(e) rep(NA_real_, n_items)
  )
  mappd_class <- vapply(mappd, .mappd_classify_scalar, character(1))

  # Effect size thresholds
  es_u_threshold  <- 0.10 * sqrt(df_uniform / 2)
  nu_es_threshold <- if (nonuniform_es == "MAPPD") 0.05 else
                       0.10 * sqrt(df_nonuniform / 2)

  # Raw omnibus p-value (used in non-uniform supplement; not BH-corrected)
  p_omnibus_raw <- stats::pchisq(chi_omnibus, df_omnibus, lower.tail = FALSE)

  # Non-uniform effect size check
  nu_es_ok <- if (nonuniform_es == "MAPPD") {
    !is.na(mappd) & mappd >= nu_es_threshold
  } else {
    !is.na(std_chi_nonuniform) & std_chi_nonuniform >= nu_es_threshold
  }

  # Primary uniform flag: BH-adjusted uniform p + effect size
  primary_uniform <- !is.na(p_uniform_adj) & p_uniform_adj < alpha &
    !is.na(std_chi_uniform) & std_chi_uniform >= es_u_threshold

  # Non-uniform supplement: bypass BH using raw omnibus p + MAPPD effect size.
  # Mirrors the LR nu_supplement approach: raw omnibus < alpha guards against
  # BH over-correction when few non-uniform items are present.  A lenient
  # gate on the marginal non-uniform p (< 2*alpha) further requires that the
  # chi_nonuniform component is at least nominally non-trivial, preventing
  # false positives driven by large chi_uniform artefacts in null items.
  # Only fires when the primary uniform criterion does not.
  nu_supplement <- !is.na(p_omnibus_raw) & p_omnibus_raw < alpha &
    !is.na(p_nonuniform) & p_nonuniform < 2 * alpha &
    nu_es_ok & !primary_uniform

  # DIF type classification (initial, before MAPPD reclassification)
  # "Uniform and Non-uniform": primary uniform flag AND marginal a-test (raw) is significant
  # "Uniform":                 primary uniform flag only
  # "Non-uniform":             nu_supplement only
  # "None":                    neither
  dif_type <- dplyr::case_when(
    is.na(p_uniform_adj)                                                  ~ NA_character_,
    primary_uniform & !is.na(p_nonuniform) & p_nonuniform < alpha         ~ "Uniform and Non-uniform",
    primary_uniform                                                        ~ "Uniform",
    nu_supplement                                                          ~ "Non-uniform",
    TRUE                                                                   ~ "None"
  )

  # Reclassify mixed items where the non-uniform ES is negligible to Uniform:
  # a significant interaction p-value with negligible MAPPD (or std_chi_nu)
  # means the non-uniform component is not practically meaningful.
  dif_type <- ifelse(
    dif_type == "Uniform and Non-uniform" & !nu_es_ok,
    "Uniform",
    dif_type
  )

  # Final type-specific flagging
  flagged_vec <- dplyr::case_when(
    dif_type == "Uniform" ~
      !is.na(p_uniform_adj) & p_uniform_adj < alpha &
      !is.na(std_chi_uniform) & std_chi_uniform >= es_u_threshold,
    dif_type == "Non-uniform" ~
      !is.na(p_nonuniform_adj) & p_nonuniform_adj < alpha &
      nu_es_ok,
    dif_type == "Uniform and Non-uniform" ~
      !is.na(p_uniform_adj) & p_uniform_adj < alpha &
      !is.na(std_chi_uniform) & std_chi_uniform >= es_u_threshold &
      nu_es_ok,
    TRUE ~ FALSE
  )

  # Type-specific primary p-values for display
  p_overall <- dplyr::case_when(
    is.na(dif_type)            ~ p_omnibus_raw,
    dif_type == "Non-uniform"  ~ p_omnibus_raw,
    TRUE                       ~ p_uniform
  )
  p_adj_primary <- dplyr::case_when(
    is.na(dif_type)            ~ p_omnibus_raw,   # raw omnibus shown for non-uniform
    dif_type == "Non-uniform"  ~ p_omnibus_raw,
    TRUE                       ~ p_uniform_adj
  )

  # ---------------------------------------------------------------------------
  # Build result data frame
  # ---------------------------------------------------------------------------

  result_df <- data.frame(
    item                 = item_cols,
    method               = "LRT",
    # Omnibus (backwards compatible)
    chi_sq               = round(chi_omnibus,    3),
    df                   = df_omnibus,
    std_chi              = round(std_chi_omnibus, 4),
    es_class             = es_class_omnibus,
    # Uniform component
    chi_uniform          = round(chi_uniform,    3),
    df_uniform           = df_uniform,
    p_uniform            = p_uniform,
    p_uniform_adj        = p_uniform_adj,
    std_chi_uniform      = round(std_chi_uniform, 4),
    es_class_uniform     = es_class_uniform,
    # Non-uniform component
    chi_nonuniform       = round(chi_nonuniform,    3),
    df_nonuniform        = df_nonuniform,
    p_nonuniform         = p_nonuniform,
    p_nonuniform_adj     = p_nonuniform_adj,
    std_chi_nonuniform   = round(std_chi_nonuniform, 4),
    # MAPPD
    mappd                = round(mappd, 4),
    mappd_class          = mappd_class,
    # Classification
    dif_type             = dif_type,
    p_overall            = p_overall,
    stringsAsFactors     = FALSE
  )

  result_df$p_adj   <- p_adj_primary
  result_df$flagged <- flagged_vec

  # ---------------------------------------------------------------------------
  # Verbose item-level output
  # ---------------------------------------------------------------------------

  if (verbose) {
    for (i in seq_len(nrow(result_df))) {
      r <- result_df[i, ]
      if (!is.na(r$flagged) && r$flagged) {
        es_part <- if (!is.na(r$dif_type) && r$dif_type == "Non-uniform") {
          if (nonuniform_es == "MAPPD")
            sprintf("MAPPD=%-8.4f", r$mappd)
          else
            sprintf("std_chi(nu)=%-6.4f", r$std_chi_nonuniform)
        } else {
          sprintf("chi2=%-8.3f  std_chi=%-6.4f", r$chi_uniform, r$std_chi_uniform)
        }
        dif_str <- if (!is.na(r$dif_type)) r$dif_type else ""
        cat(sprintf(
          "  %-20s  ⚠  %s  [%s]  %s\n",
          r$item, es_part,
          if (!is.na(r$dif_type) && r$dif_type == "Non-uniform") r$mappd_class
          else r$es_class_uniform,
          dif_str
        ))
      } else {
        cat(sprintf(
          "  %-20s  ✓  No DIF  (chi2=%.3f  p=%.3f)\n",
          r$item, r$chi_sq, r$p_adj
        ))
      }
    }
    cat("\n")
  }

  # ---------------------------------------------------------------------------
  # Direction tables for flagged items
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
      item_type <- result_df$dif_type[j]

      row_list <- list()

      # --- Difficulty direction table (Uniform and Uniform+Non-uniform) -------
      if (!is.na(item_type) &&
          item_type %in% c("Uniform", "Uniform and Non-uniform")) {

        betas    <- lls$group_betas[, j]
        ref_beta <- mean(betas, na.rm = TRUE)

        diff_rows <- lapply(seq_along(group_levels), function(gi) {
          g   <- group_levels[gi]
          dev <- if (!is.na(betas[gi])) betas[gi] - ref_beta else NA_real_
          data.frame(
            item      = item_name,
            dif_type  = "Uniform",
            group     = g,
            metric    = "Group item difficulty vs cross-group mean",
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
        row_list <- c(row_list, diff_rows)
      }

      # --- Discrimination direction table (Non-uniform and Uniform+Non-uniform)
      if (!is.na(item_type) &&
          item_type %in% c("Non-uniform", "Uniform and Non-uniform")) {

        alphas   <- lls$group_alphas[, j]
        mean_a   <- mean(alphas, na.rm = TRUE)

        disc_rows <- lapply(seq_along(group_levels), function(gi) {
          g   <- group_levels[gi]
          dev <- if (!is.na(alphas[gi])) alphas[gi] - mean_a else NA_real_
          data.frame(
            item      = item_name,
            dif_type  = "Non-uniform",
            group     = g,
            metric    = "Discrimination vs cross-group mean",
            value     = round(alphas[gi], 3),
            baseline  = round(mean_a, 3),
            deviation = round(dev, 3),
            direction = dplyr::case_when(
              is.na(dev)   ~ NA_character_,
              dev >  0.10  ~ "More discriminating than average",
              dev < -0.10  ~ "Less discriminating than average",
              dev >  0     ~ "Slightly more discriminating",
              dev <  0     ~ "Slightly less discriminating",
              TRUE         ~ "No difference"
            ),
            stringsAsFactors = FALSE
          )
        })
        row_list <- c(row_list, disc_rows)
      }

      if (length(row_list) > 0) do.call(rbind, row_list) else NULL
    })

    dir_rows_clean <- Filter(Negate(is.null), dir_rows)
    if (length(dir_rows_clean) > 0) {
      group_direction <- do.call(rbind, dir_rows_clean)
      rownames(group_direction) <- NULL
    }
  }

  .lrt_log("LRT analysis complete", verbose, success = TRUE)
  if (verbose) cat("\n")

  list(
    item_results    = result_df,
    group_direction = group_direction
  )
}


# --- Per-item constrained, intermediate, and free log-likelihoods -------------
#
# Three nested models per item, all using the anchor-based posterior:
#
#   Constrained:  a and b equal across groups      (DIF null)
#   Alpha-fixed:  a equal across groups, b free    (uniform DIF only)
#   Free:         a and b free per group            (omnibus)
#
# Returns a list with:
#   constrained  -- constrained LL per item
#   alpha_fixed  -- intermediate LL per item (a equal, b free per group)
#   free         -- free LL per item
#   group_betas  -- free-model per-group b estimates (n_groups x n_items)
#   group_alphas -- free-model per-group a estimates (n_groups x n_items)

.lrt_item_lls <- function(model, resp, group_levels, n_groups, cores = 1L) {

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

  # Worker function: compute three LLs and group params for item j
  .item_worker <- function(j) {

    x_j     <- resp_s[, j]
    obs_j   <- obs_mat[, j]
    obs_idx <- which(obs_j)

    empty <- list(ll_c_j = 0, ll_alpha_j = 0, ll_f_j = 0,
                  beta_g_j  = rep(NA_real_, n_groups),
                  alpha_g_j = rep(NA_real_, n_groups))

    if (length(obs_idx) < 5L) return(empty)

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

    eta_c   <- a_c * (nodes_obs - b_c)
    P_c     <- pmin(pmax(1 / (1 + exp(-eta_c)), 1e-10), 1 - 1e-10)
    px_c    <- ifelse(matrix(x_obs, length(obs_idx), n_nodes) == 1L, P_c, 1 - P_c)
    ll_c_j  <- sum(log(pmax(rowSums(post_obs * px_c), 1e-300)))

    # --- Intermediate: a_j shared (= a_c), b_jg free per group ---------------
    ll_j_alpha <- 0
    for (gi in seq_len(n_groups)) {
      g_idx <- which(g_obs == gi)
      if (length(g_idx) == 0L) next

      if (length(g_idx) >= 2L) {
        post_gi  <- post_obs[g_idx, , drop = FALSE]
        x_gi     <- x_obs[g_idx]
        nodes_gi <- nodes_obs[g_idx, , drop = FALSE]
        nr_ab    <- nr_update_pooled(as.double(x_gi), post_gi, nodes_gi,
                                     a_c, b_c, fix_a = TRUE, fix_b = FALSE)
        b_ag <- .clamp(nr_ab[2], -6.0, 6.0)
      } else {
        b_ag <- b_c
      }

      post_gi_all  <- post_obs[g_idx, , drop = FALSE]
      x_gi_all     <- x_obs[g_idx]
      nodes_gi_all <- nodes_obs[g_idx, , drop = FALSE]

      eta_a <- a_c * (nodes_gi_all - b_ag)
      P_a   <- pmin(pmax(1 / (1 + exp(-eta_a)), 1e-10), 1 - 1e-10)
      px_a  <- ifelse(matrix(x_gi_all, length(g_idx), n_nodes) == 1L,
                      P_a, 1 - P_a)
      ll_j_alpha <- ll_j_alpha +
        sum(log(pmax(rowSums(post_gi_all * px_a), 1e-300)))
    }

    # --- Free: a_jg, b_jg per group ------------------------------------------
    ll_j_f    <- 0
    beta_g_j  <- rep(NA_real_, n_groups)
    alpha_g_j <- rep(NA_real_, n_groups)

    for (gi in seq_len(n_groups)) {
      g_idx <- which(g_obs == gi)

      if (length(g_idx) < 2L) {
        beta_g_j[gi]  <- b_c
        alpha_g_j[gi] <- a_c
        next
      }

      post_gi  <- post_obs[g_idx, , drop = FALSE]
      x_gi     <- x_obs[g_idx]
      nodes_gi <- nodes_obs[g_idx, , drop = FALSE]

      nr_f <- nr_update_pooled(as.double(x_gi), post_gi, nodes_gi, a0, b0)
      a_f  <- .clamp(nr_f[1], 0.1, 5.0)
      b_f  <- .clamp(nr_f[2], -6.0, 6.0)
      beta_g_j[gi]  <- b_f
      alpha_g_j[gi] <- a_f

      eta_f <- a_f * (nodes_gi - b_f)
      P_f   <- pmin(pmax(1 / (1 + exp(-eta_f)), 1e-10), 1 - 1e-10)
      px_f  <- ifelse(matrix(x_gi, length(g_idx), n_nodes) == 1L, P_f, 1 - P_f)
      ll_j_f <- ll_j_f + sum(log(pmax(rowSums(post_gi * px_f), 1e-300)))
    }

    list(ll_c_j    = ll_c_j,
         ll_alpha_j = ll_j_alpha,
         ll_f_j     = ll_j_f,
         beta_g_j   = beta_g_j,
         alpha_g_j  = alpha_g_j)
  }

  # ---------------------------------------------------------------------------
  # Dispatch: serial (cores == 1) or parallel
  # ---------------------------------------------------------------------------

  item_seq <- seq_len(n_items)

  if (cores <= 1L) {

    results <- lapply(item_seq, .item_worker)

  } else if (.Platform$OS.type == "windows") {

    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Load the package (dev version) on each worker
    pkg_path <- system.file(package = "iDIFr")
    if (nchar(pkg_path) == 0L) {
      # Under devtools::load_all() the package is not installed — load from source
      src_path <- tryCatch(
        find.package("iDIFr"),
        error = function(e) NULL
      )
      if (is.null(src_path)) {
        # Fall back to the working directory heuristic
        src_path <- getwd()
      }
      parallel::clusterCall(cl, function(p) devtools::load_all(p), src_path)
    } else {
      parallel::clusterEvalQ(cl, library(iDIFr))
    }

    # Export the environment objects the worker closure captures
    parallel::clusterExport(
      cl,
      varlist = c("resp_s", "obs_mat", "post", "nodes_i", "g_int",
                  "n_groups", "n_nodes"),
      envir   = environment()
    )

    results <- parallel::parLapply(cl, item_seq, .item_worker)

  } else {

    results <- parallel::mclapply(item_seq, .item_worker, mc.cores = cores)

  }

  # ---------------------------------------------------------------------------
  # Unpack results list into arrays
  # ---------------------------------------------------------------------------

  ll_c         <- numeric(n_items)
  ll_alpha     <- numeric(n_items)
  ll_f         <- numeric(n_items)
  betas_g_mat  <- matrix(NA_real_, nrow = n_groups, ncol = n_items)
  alphas_g_mat <- matrix(NA_real_, nrow = n_groups, ncol = n_items)

  for (j in item_seq) {
    r               <- results[[j]]
    ll_c[j]         <- r$ll_c_j
    ll_alpha[j]     <- r$ll_alpha_j
    ll_f[j]         <- r$ll_f_j
    betas_g_mat[, j]  <- r$beta_g_j
    alphas_g_mat[, j] <- r$alpha_g_j
  }

  list(
    constrained  = ll_c,
    alpha_fixed  = ll_alpha,
    free         = ll_f,
    group_betas  = betas_g_mat,
    group_alphas = alphas_g_mat
  )
}


# --- MAPPD for IRT models -------------------------------------------------------
#
# Maximum Absolute Predicted Probability Difference for the free IRT model.
# a_mat: n_groups x n_items matrix of discrimination parameters
# b_mat: n_groups x n_items matrix of difficulty parameters
# Returns a numeric vector of length n_items.

.mappd_irt <- function(a_mat, b_mat, group_levels, n_grid = 100) {
  theta_grid <- seq(-3, 3, length.out = n_grid)
  n_groups   <- nrow(a_mat)
  n_items    <- ncol(a_mat)

  # Vectorised over grid points and group pairs using matrix operations.
  # For item j: build a (n_grid x n_groups) probability matrix, then find the
  # max absolute pairwise difference across groups in one pass.
  mappd_vec <- vapply(seq_len(n_items), function(j) {

    if (anyNA(a_mat[, j]) || anyNA(b_mat[, j])) return(NA_real_)

    # P[t, g] = 1 / (1 + exp(-a_g * (theta_t - b_g)))
    # theta_grid is n_grid; a_mat[,j] and b_mat[,j] are n_groups vectors.
    # Use outer to build the n_grid x n_groups eta matrix.
    a_j   <- a_mat[, j]   # length n_groups
    b_j   <- b_mat[, j]   # length n_groups
    # eta[t,g] = a_j[g] * (theta_grid[t] - b_j[g])
    eta   <- outer(theta_grid, a_j) - matrix(a_j * b_j, nrow = n_grid,
                                             ncol = n_groups, byrow = TRUE)
    p_mat <- 1 / (1 + exp(-eta))   # n_grid x n_groups

    # All pairs (i, k) with i < k: vectorise via column subtraction
    max_diff <- 0
    for (i in seq_len(n_groups - 1L)) {
      for (k in seq.int(i + 1L, n_groups)) {
        d <- max(abs(p_mat[, i] - p_mat[, k]))
        if (d > max_diff) max_diff <- d
      }
    }
    max_diff

  }, numeric(1))

  # --- Verification on construction (only when run interactively / testing) ---
  # Uncomment the block below to cross-check against the old loop-based result.
  #
  # old_result <- .mappd_irt_loop(a_mat, b_mat, group_levels, n_grid)
  # stopifnot(all.equal(mappd_vec, old_result))

  mappd_vec
}

# Self-contained verification: run once on a tiny known case and print result.
local({
  a_test <- matrix(c(1, 1.5, 0.8, 1.2), nrow = 2, ncol = 2)  # 2 groups x 2 items
  b_test <- matrix(c(0, 0.5, -0.5, 0.3), nrow = 2, ncol = 2)
  gl     <- c("G1", "G2")

  # Reference: brute-force loop
  ref_loop <- function(a_mat, b_mat, n_grid = 100) {
    tg <- seq(-3, 3, length.out = n_grid)
    ng <- nrow(a_mat); ni <- ncol(a_mat)
    out <- numeric(ni)
    for (j in seq_len(ni)) {
      if (anyNA(a_mat[, j]) || anyNA(b_mat[, j])) { out[j] <- NA_real_; next }
      pm <- matrix(0, n_grid, ng)
      for (gi in seq_len(ng)) pm[, gi] <- 1/(1+exp(-a_mat[gi,j]*(tg - b_mat[gi,j])))
      mx <- 0
      for (i in seq_len(ng)) for (k in seq_len(ng)) {
        if (i >= k) next
        mx <- max(mx, max(abs(pm[,i]-pm[,k])))
      }
      out[j] <- mx
    }
    out
  }

  vec_result  <- .mappd_irt(a_test, b_test, gl)
  loop_result <- ref_loop(a_test, b_test)
  ok <- isTRUE(all.equal(vec_result, loop_result))
  cat(sprintf("[.mappd_irt vectorisation check] identical to loop: %s\n", ok))
})

# Scalar wrapper for vapply compatibility
.mappd_classify_scalar <- function(x) {
  if (is.na(x))    return(NA_character_)
  if (x >= 0.10)   return("Large")
  if (x >= 0.05)   return("Moderate")
  "Negligible"
}


# --- Logging helper -----------------------------------------------------------

.lrt_log <- function(msg, verbose, detail = NULL, success = FALSE) {
  if (!verbose) return(invisible(NULL))
  timestamp <- format(Sys.time(), "[%H:%M:%S]")
  symbol    <- if (success) "✓" else " "
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
    item                 = item_name,
    method               = "LRT",
    chi_sq               = NA_real_,
    df                   = NA_integer_,
    std_chi              = NA_real_,
    es_class             = NA_character_,
    chi_uniform          = NA_real_,
    df_uniform           = NA_integer_,
    p_uniform            = NA_real_,
    p_uniform_adj        = NA_real_,
    std_chi_uniform      = NA_real_,
    es_class_uniform     = NA_character_,
    chi_nonuniform       = NA_real_,
    df_nonuniform        = NA_integer_,
    p_nonuniform         = NA_real_,
    p_nonuniform_adj     = NA_real_,
    std_chi_nonuniform   = NA_real_,
    mappd                = NA_real_,
    mappd_class          = NA_character_,
    dif_type             = NA_character_,
    p_overall            = NA_real_,
    p_adj                = NA_real_,
    flagged              = NA,
    stringsAsFactors     = FALSE
  )
}

.lrt_na_results <- function(item_cols) {
  list(
    item_results    = do.call(rbind, lapply(item_cols, .lrt_na_row)),
    group_direction = NULL
  )
}

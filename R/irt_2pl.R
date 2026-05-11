# =============================================================================
# irt_2pl.R
# Marginal Maximum Likelihood 2PL IRT estimation via EM algorithm.
#
# E-step and pooled NR M-step are implemented in C++ (src/em_2pl.cpp).
# Single group: < 0.5s for n=1000, 20 items.
# Multigroup:   < 30s for n=1000, 20 items, 12 groups.
#
# Parameterisation: P(X=1|theta,a,b) = 1/(1 + exp(-a*(theta-b)))
# Quadrature: Gauss-Hermite, 21 nodes, scaled per group N(mu, sigma2)
# Reference:  Bock & Aitkin (1981) Psychometrika 46:443-459
# =============================================================================

#' Fit a 2PL IRT model via marginal maximum likelihood (EM)
#'
#' @param resp      Integer matrix (0/1/NA). Rows=persons, cols=items.
#' @param group     Character/factor vector of group membership (length=nrow(resp)).
#'                  NULL for single-group calibration.
#' @param constrain Parameter constraint across groups:
#'   \itemize{
#'     \item \code{"items"} — a and b equal across groups (DIF null hypothesis)
#'     \item \code{"none"}  — all parameters free (separate group calibrations)
#'     \item \code{"alpha"} — a fixed across groups, b free (uniform DIF test)
#'     \item \code{"beta"}  — b fixed across groups, a free (non-uniform DIF test)
#'   }
#' @param n_nodes   Number of quadrature nodes. Default 21.
#' @param max_iter  Maximum EM iterations. Default 500.
#' @param tol       Convergence tolerance on log-likelihood change. Default 1e-4.
#' @param verbose   Print iteration log. Default FALSE.
#'
#' @return Object of class \code{irt_2pl}.
#' @export
fit_2pl <- function(resp,
                    group     = NULL,
                    constrain = "items",
                    n_nodes   = 21,
                    max_iter  = 500,
                    tol       = 1e-4,
                    verbose   = FALSE) {

  resp <- as.matrix(resp)
  stopifnot(all(resp %in% c(0L, 1L, NA)))

  valid_c <- c("items", "none", "alpha", "beta")
  if (!constrain %in% valid_c)
    stop("constrain must be one of: ", paste(valid_c, collapse=", "), call.=FALSE)

  n_persons  <- nrow(resp)
  n_items    <- ncol(resp)
  item_names <- if (!is.null(colnames(resp))) colnames(resp) else
                paste0("item_", seq_len(n_items))
  colnames(resp) <- item_names

  # NA indicator matrix — TRUE where observed
  obs_mat <- !is.na(resp)
  resp_safe <- resp; resp_safe[!obs_mat] <- 0L  # fill NA with 0 for maths

  # --- Groups ----------------------------------------------------------------

  if (is.null(group)) group <- rep("G1", n_persons)
  group        <- as.character(group)
  group_levels <- sort(unique(group))
  n_groups     <- length(group_levels)
  g_int        <- match(group, group_levels)

  # --- Quadrature ------------------------------------------------------------

  gh          <- .gh_nodes(n_nodes)
  std_nodes   <- gh$nodes    # length n_nodes, N(0,1) scale
  std_weights <- gh$weights  # sum to 1

  # nodes matrix: n_nodes x 1 (broadcast across persons later)
  # theta_gk = mu_g + sqrt(sigma2_g) * std_nodes[k]

  # --- Initialise item parameters --------------------------------------------

  # a in (0.1, 5), b in (-6, 6)
  # Initial b from empirical p-correct, a=1
  pc <- colMeans(resp_safe * obs_mat, na.rm=TRUE) /
        pmax(colMeans(obs_mat), 0.01)
  pc <- pmin(pmax(pc, 0.01), 0.99)

  # a_mat[g, j], b_mat[g, j]
  a_mat <- matrix(1.0,          nrow=n_groups, ncol=n_items,
                  dimnames=list(group_levels, item_names))
  b_mat <- matrix(qlogis(1-pc), nrow=n_groups, ncol=n_items,
                  byrow=TRUE,   dimnames=list(group_levels, item_names))

  # Group ability: mu (free for g>1), sigma2 (free for all)
  # Reference group: mu=0, sigma2=1 (identification)
  mu_g     <- c(0, rep(0, n_groups-1)); names(mu_g)     <- group_levels
  sigma2_g <- rep(1, n_groups);         names(sigma2_g) <- group_levels

  # --- EM loop ---------------------------------------------------------------

  ll_hist   <- numeric(max_iter)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    # ---- E step -----------------------------------------------------------

    # Compute scaled nodes per group: theta_gk = mu_g + sqrt(s2_g)*z_k
    # nodes_g: n_groups x n_nodes
    nodes_g <- outer(sqrt(sigma2_g), std_nodes, "*") +
               matrix(mu_g, nrow=n_groups, ncol=n_nodes)

    # Per-person parameter matrices (constant across nodes)
    a_i <- a_mat[g_int, , drop = FALSE]   # n_persons x n_items
    b_i <- b_mat[g_int, , drop = FALSE]

    # C++ kernel: builds n_persons x n_nodes log-lik matrix
    loglik_mat <- compute_log_lik_matrix(
      resp_safe, obs_mat, a_i, b_i, nodes_g, g_int
    )

    # Posterior: P(theta_k | x_i) proportional to weight_k * exp(loglik_mat[i,k])
    log_post <- sweep(loglik_mat, 2, log(std_weights), "+")
    log_norm <- .lse_rows(log_post)                    # n_persons log-normaliser
    loglik   <- sum(log_norm)
    post     <- exp(log_post - log_norm)               # n_persons x n_nodes

    ll_hist[iter] <- loglik

    if (verbose) {
      delta_str <- if (iter > 1) sprintf("  delta=%.6f", loglik - ll_hist[iter-1]) else ""
      cat(sprintf("  iter %3d  loglik=%.4f%s\n", iter, loglik, delta_str))
    }

    if (iter > 1 && abs(loglik - ll_hist[iter-1]) < tol) {
      converged <- TRUE; break
    }

    # ---- M step -----------------------------------------------------------

    # Sufficient statistics per group per node:
    # r_gjk = sum_{i in g} post[i,k] * x_ij   (expected correct at node k)
    # f_gjk = sum_{i in g} post[i,k]            (expected count at node k)
    #
    # For constrained params (shared across groups):
    # r_jk  = sum_g r_gjk,  f_jk = sum_g f_gjk  (pool across groups)

    new_a <- a_mat
    new_b <- b_mat

    # Update item parameters
    for (j in seq_len(n_items)) {

      x_j   <- resp_safe[, j]
      obs_j <- obs_mat[, j]

      if (constrain == "none") {
        # Separate update per group
        for (gi in seq_len(n_groups)) {
          idx    <- g_int == gi & obs_j
          if (sum(idx) == 0) next
          post_g <- post[idx, , drop=FALSE]
          x_g    <- x_j[idx]
          nd     <- nodes_g[gi, ]

          nr <- .nr_2pl(x_g, post_g, nd,
                        a_mat[gi,j], b_mat[gi,j],
                        fix_a=FALSE, fix_b=FALSE)
          new_a[gi,j] <- .clamp(nr$a, 0.1, 5.0)
          new_b[gi,j] <- .clamp(nr$b, -6.0, 6.0)
        }

      } else {
        # Pooled update using all persons
        # Build per-person node vectors (their group's scaled nodes at each k)
        # nodes_ik: n_persons_obs x n_nodes
        obs_idx  <- which(obs_j)
        post_obs <- post[obs_idx, , drop=FALSE]
        x_obs    <- x_j[obs_idx]
        g_obs    <- g_int[obs_idx]

        # Scaled nodes for each observed person
        nodes_ik <- nodes_g[g_obs, , drop=FALSE]  # n_obs x n_nodes

        fix_a <- constrain == "beta"
        fix_b <- constrain == "alpha"

        nr <- nr_update_pooled(as.double(x_obs), post_obs, nodes_ik,
                               a_mat[1, j], b_mat[1, j],
                               fix_a = fix_a, fix_b = fix_b)

        new_a_j <- .clamp(nr[1], 0.1, 5.0)
        new_b_j <- .clamp(nr[2], -6.0, 6.0)

        for (gi in seq_len(n_groups)) {
          if (!fix_a) new_a[gi,j] <- new_a_j
          if (!fix_b) new_b[gi,j] <- new_b_j
        }
      }
    }

    a_mat <- new_a
    b_mat <- new_b

    # Update group ability distributions (free groups only, g > 1)
    for (gi in seq_len(n_groups)) {
      if (gi == 1) next
      idx     <- which(g_int == gi)
      post_g  <- post[idx, , drop=FALSE]
      nd      <- nodes_g[gi, ]  # n_nodes

      # E[theta] and E[theta^2] under posterior
      E_t  <- rowSums(post_g * matrix(nd,   nrow=length(idx), ncol=n_nodes, byrow=TRUE))
      E_t2 <- rowSums(post_g * matrix(nd^2, nrow=length(idx), ncol=n_nodes, byrow=TRUE))

      mu_g[gi]     <- mean(E_t)
      sigma2_g[gi] <- max(0.01, mean(E_t2) - mean(E_t)^2)
    }

  }  # end EM

  if (!converged && verbose)
    warning(sprintf("EM did not converge in %d iterations.", max_iter))

  # --- Output ----------------------------------------------------------------

  if (constrain %in% c("items", "alpha", "beta")) {
    item_params <- data.frame(item=item_names,
                               a=a_mat[1,], b=b_mat[1,],
                               stringsAsFactors=FALSE, row.names=NULL)
  } else {
    item_params <- data.frame(item=item_names,
                               a=a_mat[1,], b=b_mat[1,],
                               stringsAsFactors=FALSE, row.names=NULL)
    for (gi in seq_len(n_groups)) {
      g <- group_levels[gi]
      item_params[[paste0("a_",g)]] <- a_mat[gi,]
      item_params[[paste0("b_",g)]] <- b_mat[gi,]
    }
  }

  group_params <- data.frame(group=group_levels, mu=mu_g, sigma2=sigma2_g,
                              stringsAsFactors=FALSE, row.names=NULL)

  structure(
    list(
      item_params    = item_params,
      group_params   = group_params,
      loglik         = ll_hist[iter],
      loglik_history = ll_hist[seq_len(iter)],
      iter           = iter,
      converged      = converged,
      n_items        = n_items,
      n_persons      = n_persons,
      n_groups       = n_groups,
      group_levels   = group_levels,
      constrain      = constrain,
      nodes          = std_nodes,
      weights        = std_weights,
      posterior      = post,
      a_mat          = a_mat,
      b_mat          = b_mat,
      mu_vec         = mu_g,
      sigma2_vec     = sigma2_g,
      resp           = resp,
      group_vector   = group,
      g_int          = g_int
    ),
    class = "irt_2pl"
  )
}


# --- S3 methods ---------------------------------------------------------------

#' @export
print.irt_2pl <- function(x, ...) {
  cat("2PL IRT Model\n")
  cat(sprintf("  Items:      %d\n", x$n_items))
  cat(sprintf("  Persons:    %d\n", x$n_persons))
  cat(sprintf("  Groups:     %d\n", x$n_groups))
  cat(sprintf("  Constraint: %s\n", x$constrain))
  cat(sprintf("  Log-lik:    %.4f\n", x$loglik))
  cat(sprintf("  Converged:  %s  (%d iterations)\n\n", x$converged, x$iter))
  cat("Item parameters:\n")
  print(round(x$item_params, 4))
  if (x$n_groups > 1) {
    cat("\nGroup ability parameters:\n")
    print(round(x$group_params, 4))
  }
  invisible(x)
}

#' @export
logLik.irt_2pl <- function(object, ...) {
  structure(object$loglik, df=NA_integer_, class="logLik")
}


# --- Per-item log-likelihood decomposition ------------------------------------

#' Compute per-item log-likelihood contributions from a fitted irt_2pl model
#'
#' Uses local independence to decompose total LL into item contributions:
#' \code{LL_j = sum_i log P(x_ij | posterior_i)}
#' where \code{P(x_ij | posterior_i) = sum_k posterior[i,k] * P(x_ij | theta_k)}
#'
#' @param model    An \code{irt_2pl} object.
#' @param resp     Response matrix (0/1/NA). Defaults to model$resp.
#' @param post     Posterior matrix (persons x nodes). Defaults to model$posterior.
#' @param gi       Group index (integer). Used to select group-specific item
#'                 params and ability nodes. Default 1.
#'
#' @return Numeric vector of length n_items.
#' @export
item_loglik <- function(model, resp=NULL, post=NULL, gi=1) {

  if (is.null(resp)) resp <- model$resp
  if (is.null(post)) post <- model$posterior

  resp     <- as.matrix(resp)
  n_p      <- nrow(resp)
  n_items  <- ncol(resp)
  n_nodes  <- ncol(post)

  a_j <- model$a_mat[gi, ]
  b_j <- model$b_mat[gi, ]
  mu  <- model$mu_vec[gi]
  s2  <- model$sigma2_vec[gi]
  nd  <- mu + sqrt(s2) * model$nodes   # scaled nodes for this group

  # P matrix: n_items x n_nodes
  # P[j,k] = P(X=1 | theta_k, a_j, b_j)
  eta <- outer(a_j, nd) - a_j * b_j   # n_items x n_nodes  (a*(theta-b))
  P   <- 1 / (1 + exp(-eta))
  P   <- pmin(pmax(P, 1e-10), 1-1e-10)

  item_ll <- numeric(n_items)
  obs_mat <- !is.na(resp)
  resp_s  <- resp; resp_s[!obs_mat] <- 0L

  for (j in seq_len(n_items)) {
    # For each person: P(x_ij | post_i) = sum_k post[i,k] * P[j,k]^x * (1-P[j,k])^(1-x)
    # Vectorised over persons:
    p_jk     <- P[j, ]                    # n_nodes
    x_j      <- resp_s[, j]               # n_persons
    obs_j    <- obs_mat[, j]

    # P(x_ij | theta_k): n_persons x n_nodes
    px <- matrix(p_jk, nrow=n_p, ncol=n_nodes, byrow=TRUE)
    px <- ifelse(matrix(x_j, n_p, n_nodes) == 1, px, 1-px)
    px[!obs_j, ] <- 1  # unobserved contributes 1 (no information)

    # Integrate over posterior
    px_post  <- rowSums(post * px)
    item_ll[j] <- sum(log(pmax(px_post, 1e-300))[obs_j])
  }

  item_ll
}


#' Per-item LL for a multigroup constrained model
#'
#' For the constrained model, each person uses shared item params but
#' their own group-specific ability nodes.
#'
#' @param model   An \code{irt_2pl} object with constrain != "none".
#' @param resp    Response matrix. Defaults to model$resp.
#' @param post    Posterior matrix. Defaults to model$posterior.
#' @return Numeric vector of length n_items.
#' @export
item_loglik_mg <- function(model, resp=NULL, post=NULL) {

  if (is.null(resp)) resp <- model$resp
  if (is.null(post)) post <- model$posterior

  resp     <- as.matrix(resp)
  n_p      <- nrow(resp)
  n_items  <- ncol(resp)
  n_nodes  <- ncol(post)
  n_groups <- model$n_groups
  g_int    <- model$g_int

  # Shared item params (row 1)
  a_j <- model$a_mat[1, ]
  b_j <- model$b_mat[1, ]

  # Group-scaled nodes: n_groups x n_nodes
  nodes_g <- outer(sqrt(model$sigma2_vec), model$nodes, "*") +
             matrix(model$mu_vec, nrow=n_groups, ncol=n_nodes)

  # Per-person nodes: n_persons x n_nodes
  nodes_i <- nodes_g[g_int, , drop=FALSE]

  obs_mat <- !is.na(resp)
  resp_s  <- resp; resp_s[!obs_mat] <- 0L

  item_ll <- numeric(n_items)

  for (j in seq_len(n_items)) {
    # eta_ik = a_j * (theta_ik - b_j): n_persons x n_nodes
    eta_ik <- a_j[j] * (nodes_i - b_j[j])   # broadcasting
    P_ik   <- 1 / (1 + exp(-eta_ik))
    P_ik   <- pmin(pmax(P_ik, 1e-10), 1-1e-10)

    x_j  <- resp_s[, j]
    obs_j <- obs_mat[, j]

    # P(x_ij | theta_k): n_persons x n_nodes
    px <- ifelse(matrix(x_j, n_p, n_nodes) == 1, P_ik, 1-P_ik)
    px[!obs_j, ] <- 1

    px_post  <- rowSums(post * px)
    item_ll[j] <- sum(log(pmax(px_post, 1e-300))[obs_j])
  }

  item_ll
}


# --- Quadrature ---------------------------------------------------------------


.gh_nodes <- function(n) {
  b   <- sqrt(seq_len(n - 1) / 2)
  J   <- diag(0, n)
  for (k in seq_len(n - 1)) {
    J[k,   k+1] <- b[k]
    J[k+1, k  ] <- b[k]
  }
  
  eig   <- eigen(J, symmetric = TRUE)
  idx   <- order(eig$values)
  nodes <- eig$values[idx] * sqrt(2)   # scale to N(0,1): var=1 not 0.5
  
  w_raw <- eig$vectors[1, idx]^2 * sqrt(2 * pi)
  w     <- w_raw / sum(w_raw)
  
  list(nodes = nodes, weights = w)
}


# --- Newton-Raphson for single group ------------------------------------------

.nr_2pl <- function(x, post, nodes, a, b,
                     fix_a=FALSE, fix_b=FALSE,
                     max_nr=20, tol_nr=1e-5) {

  # x: n_obs vector of responses
  # post: n_obs x n_nodes posterior
  # nodes: n_nodes vector of quadrature points (group-scaled)
  n_nodes <- length(nodes)

  # Sufficient statistics
  # r_k = sum_i post[i,k] * x_i   (expected correct at node k)
  # f_k = sum_i post[i,k]          (expected total at node k)
  r_k <- colSums(post * x)
  f_k <- colSums(post)

  for (nr in seq_len(max_nr)) {

    p_k <- 1/(1+exp(-a*(nodes-b)));  p_k <- pmin(pmax(p_k,1e-10),1-1e-10)
    d_k <- p_k*(1-p_k)
    z_k <- nodes-b
    res <- r_k - f_k*p_k

    g_a <- if(!fix_a) sum(res*z_k)          else 0
    h_a <- if(!fix_a) -sum(f_k*d_k*z_k^2)  else -1
    g_b <- if(!fix_b) sum(res*(-a))         else 0
    h_b <- if(!fix_b) -sum(f_k*d_k*a^2)    else -1

    if (!fix_a && abs(h_a)>1e-10) a <- a - g_a/h_a
    if (!fix_b && abs(h_b)>1e-10) b <- b - g_b/h_b

    if (abs(g_a)+abs(g_b) < tol_nr) break
  }

  list(a=a, b=b)
}


# --- Utilities ----------------------------------------------------------------

.lse_rows <- function(mat) {
  # Numerically stable log-sum-exp across columns (per row)
  m <- apply(mat, 1, max)
  m + log(rowSums(exp(mat - m)))
}

.clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

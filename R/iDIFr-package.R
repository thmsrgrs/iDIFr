#' iDIFr: Intersectional Differential Item Functioning Analysis in R
#'
#' @description
#' A user-friendly toolkit for detecting Differential Item Functioning (DIF)
#' using Logistic Regression (LR), the IRT Likelihood Ratio Test (LRT), and
#' Intersectional Decomposition (ID). Designed for both standard two-group
#' and intersectional multi-group designs.
#'
#' ## Key functions
#'
#' - [idifr()]: Main entry point — run DIF analysis
#' - [check_groups()]: Explore group structure and cell sizes
#' - [merge_groups()]: Combine sparse intersectional cells
#' - [tidy.idifr()]: Extract results as a flat data frame
#'
#' ## Quick start
#'
#' ```r
#' library(iDIFr)
#'
#' # Check your group structure first
#' check_groups(my_data, group = ~ gender * nationality * age_band)
#'
#' # Run DIF analysis
#' result <- idifr(
#'   data   = my_data,
#'   items  = 1:20,
#'   group  = ~ gender * nationality * age_band,
#'   method = c("LR", "LRT")
#' )
#'
#' print(result)    # Flagged items with effect sizes
#' summary(result)  # Full breakdown by method
#' plot(result)     # Effect size heatmap
#' tidy(result)     # Flat data frame for further analysis
#' ```
#'
#' @keywords internal
#' @importFrom stats qlogis reorder p.adjust pchisq glm binomial anova qnorm
#' @importFrom rlang .data
#' @importFrom dplyr case_when mutate filter group_by summarise left_join arrange desc
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 theme_minimal labs theme scale_fill_manual facet_wrap element_blank element_text
#' @importFrom cli cli_alert_warning cli_alert_info cli_alert_success cli_h1 cli_h2 cli_text col_yellow col_red col_green
#' @importFrom strucchange efp sctest
#' @importFrom Rcpp sourceCpp
#' @useDynLib iDIFr, .registration = TRUE
"_PACKAGE"

utils::globalVariables(c(".data", "n", "method", "item", "es", "flagged"))


#' Generate synthetic DIF data for testing and simulation
#'
#' @description
#' A helper function for generating synthetic dichotomous item response data
#' with known DIF structure. Useful for testing the package, running pilot
#' analyses, and simulation studies.
#'
#' @param n_persons Integer. Total number of respondents.
#' @param n_items Integer. Number of items. Default is 20.
#' @param n_groups Integer. Number of groups. Default is 2.
#' @param dif_items Integer vector. Which items have DIF. Default is `c(3, 7)`.
#' @param dif_effect Numeric. Size of DIF effect in logits. Default is `0.8`.
#' @param dif_type Character. `"uniform"` (default) or `"nonuniform"`.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A data frame with item response columns (`item_1`, `item_2`, ...) and
#'   a `group` column.
#'
#' @examples
#' \dontrun{
#' # Generate a simple 2-group dataset
#' dat <- simulate_dif(n_persons = 500, n_items = 20, dif_items = c(3, 7))
#'
#' # Check it works with idifr
#' result <- idifr(dat, items = 1:20, group = ~ group, method = "LR")
#' }
#'
#' @export
simulate_dif <- function(n_persons  = 500,
                         n_items    = 20,
                         n_groups   = 2,
                         dif_items  = c(3, 7),
                         dif_effect = 0.8,
                         dif_type   = "uniform",
                         seed       = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # --- Group assignment -------------------------------------------------------

  group        <- factor(sample(paste0("G", seq_len(n_groups)),
                                n_persons, replace = TRUE))
  group_int    <- as.integer(group)   # G1=1, G2=2, ...
  is_focal     <- group_int > 1       # all non-reference groups are focal

  # --- Person abilities -------------------------------------------------------

  theta <- stats::rnorm(n_persons, mean = 0, sd = 1)

  # --- Item parameters --------------------------------------------------------
  # ONE set of parameters shared across ALL groups.
  # DIF is introduced by modifying parameters for focal groups on DIF items
  # only — non-DIF items use identical parameters for every group.

  a <- stats::runif(n_items, 0.5, 2.0)   # discrimination (shared)
  b <- stats::rnorm(n_items, 0, 1)       # difficulty     (shared)

  dif_items <- dif_items[dif_items >= 1L & dif_items <= n_items]

  # --- Simulate responses vectorised -----------------------------------------
  # Build n_persons x n_items matrices of effective a and b parameters,
  # then compute P and draw responses in one vectorised step.

  # Start with shared parameters broadcast to all persons
  a_mat <- matrix(a, nrow = n_persons, ncol = n_items, byrow = TRUE)
  b_mat <- matrix(b, nrow = n_persons, ncol = n_items, byrow = TRUE)

  # Apply DIF modifications for focal group persons on DIF items
  if (length(dif_items) > 0) {
    for (j in dif_items) {
      if (dif_type == "uniform") {
        # Shift difficulty for focal group only — discrimination unchanged
        b_mat[is_focal, j] <- b[j] + dif_effect

      } else if (dif_type == "nonuniform") {
        # Shift both discrimination and difficulty for focal group
        a_mat[is_focal, j] <- a[j] * (1 + dif_effect * 0.5)
        b_mat[is_focal, j] <- b[j] + dif_effect * 0.3
      }
    }
  }

  # Compute P(X=1 | theta, a, b) for all persons and items simultaneously
  # theta_mat: n_persons x n_items (each row is the same theta value)
  theta_mat <- matrix(theta, nrow = n_persons, ncol = n_items)
  P         <- 1 / (1 + exp(-a_mat * (theta_mat - b_mat)))

  # Draw binary responses
  responses <- matrix(
    stats::rbinom(n_persons * n_items, 1, as.vector(P)),
    nrow = n_persons, ncol = n_items
  )

  # --- Build output data frame ------------------------------------------------

  df       <- as.data.frame(responses)
  names(df) <- paste0("item_", seq_len(n_items))
  df$group  <- group

  # Store true parameters as attributes for validation use
  attr(df, "true_a")       <- a
  attr(df, "true_b")       <- b
  attr(df, "dif_items")    <- dif_items
  attr(df, "dif_effect")   <- dif_effect
  attr(df, "dif_type")     <- dif_type

  df
}

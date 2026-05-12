#' iDIFr: Intersectional Differential Item Functioning Analysis in R
#'
#' @description
#' A user-friendly toolkit for detecting Differential Item Functioning (DIF)
#' using Logistic Regression (LR), the IRT Likelihood Ratio Test (LRT), and
#' model-based recursive partitioning (MOB). Designed for both standard
#' two-group and intersectional multi-group designs, with built-in
#' Intersectional Contrast Analysis (ICA) via the `ica = TRUE` argument.
#'
#' ## Key functions
#'
#' - [idifr()]: Main entry point -- run DIF analysis (set `ica = TRUE` for ICA)
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
#' @importFrom stats qlogis reorder p.adjust pchisq glm binomial anova qnorm setNames
#' @importFrom utils combn head
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
#' Generates synthetic dichotomous item response data with a known DIF
#' structure. Supports three DIF patterns: standard group DIF (`"standard"`),
#' DIF confined to a single intersectional cell (`"intersection"`), and a
#' mixture of both (`"mixed"`).
#'
#' @param n_persons Integer. Total number of respondents.
#' @param n_items Integer. Number of items. Default 20.
#' @param n_groups Integer. Number of groups. Default 2.
#' @param dif_items Which items have DIF. For `dif_structure = "standard"` or
#'   `"intersection"`, an integer vector (e.g. `c(3, 7)`). For
#'   `dif_structure = "mixed"`, either a named list
#'   `list(standard = c(3,7), intersection = c(12,15))` or a plain integer
#'   vector (first half gets standard DIF, second half intersection DIF).
#'   Default `c(3, 7)`.
#' @param dif_effect Numeric. DIF shift size in logits. Default `0.8`.
#' @param dif_type `"uniform"` (difficulty shift only, default) or
#'   `"nonuniform"` (both difficulty and discrimination shifted).
#' @param dif_structure One of `"standard"` (default), `"intersection"`, or
#'   `"mixed"`. `"standard"` replicates the original behaviour. `"intersection"`
#'   applies DIF only to the specific intersectional cell in `dif_group`.
#'   `"mixed"` applies standard DIF to some items and intersection DIF to
#'   others.
#' @param dif_group Named list identifying the target intersectional cell for
#'   intersection DIF. Variable names must match `demo_vars` or `"group"`.
#'   Example: `list(group = "G1", nationality = "UK", age_band = "Young")`.
#'   Required when `dif_structure` is `"intersection"` or `"mixed"`.
#' @param demo_vars Named list of additional demographic variables to add, with
#'   their levels. Persons are assigned randomly with uniform probability.
#'   Example: `list(nationality = c("UK", "DE", "FR"), age_band = c("Young",
#'   "Old"))`. Required when `dif_structure` is `"intersection"` or `"mixed"`.
#' @param seed Integer random seed for reproducibility.
#'
#' @return A data frame with item response columns (`item_1`, `item_2`, ...),
#'   a `group` column, and any additional columns specified in `demo_vars`.
#'   True item parameters and DIF metadata are stored as attributes.
#'
#' @examples
#' \dontrun{
#' # Standard DIF — unchanged behaviour
#' dat <- simulate_dif(500, 20, 2, c(3, 7), 1.0)
#'
#' # Intersection-only DIF
#' dat_ix <- simulate_dif(
#'   n_persons     = 2000,
#'   n_items       = 20,
#'   dif_items     = c(5, 12),
#'   dif_effect    = 1.5,
#'   dif_structure = "intersection",
#'   dif_group     = list(group = "G1", nationality = "UK", age_band = "Young"),
#'   demo_vars     = list(nationality = c("UK", "DE", "FR"),
#'                        age_band    = c("Young", "Old")),
#'   seed          = 42
#' )
#'
#' # Mixed DIF
#' dat_mix <- simulate_dif(
#'   n_persons     = 2000,
#'   n_items       = 20,
#'   dif_items     = list(standard = c(3, 7), intersection = c(12, 15)),
#'   dif_effect    = 1.0,
#'   dif_structure = "mixed",
#'   dif_group     = list(group = "G1", nationality = "UK", age_band = "Young"),
#'   demo_vars     = list(nationality = c("UK", "DE", "FR"),
#'                        age_band    = c("Young", "Old")),
#'   seed          = 42
#' )
#' }
#'
#' @export
simulate_dif <- function(n_persons     = 500,
                         n_items       = 20,
                         n_groups      = 2,
                         dif_items     = c(3, 7),
                         dif_effect    = 0.8,
                         dif_type      = "uniform",
                         dif_structure = "standard",
                         dif_group     = NULL,
                         demo_vars     = NULL,
                         seed          = NULL) {

  dif_structure <- match.arg(dif_structure, c("standard", "intersection", "mixed"))
  dif_type      <- match.arg(dif_type,      c("uniform",  "nonuniform"))

  if (!is.null(seed)) set.seed(seed)

  # --- Validate intersection arguments ----------------------------------------

  if (dif_structure %in% c("intersection", "mixed")) {
    if (is.null(demo_vars)) {
      stop("`demo_vars` must be specified when dif_structure = '",
           dif_structure, "'.", call. = FALSE)
    }
    if (is.null(dif_group)) {
      stop("`dif_group` must be specified when dif_structure = '",
           dif_structure, "'.", call. = FALSE)
    }
    all_demo_names <- c("group", names(demo_vars))
    bad_vars <- setdiff(names(dif_group), all_demo_names)
    if (length(bad_vars) > 0) {
      stop("`dif_group` references unknown variable(s): ",
           paste(bad_vars, collapse = ", "), call. = FALSE)
    }
  }

  # --- Parse dif_items --------------------------------------------------------

  if (dif_structure == "mixed") {
    if (is.list(dif_items)) {
      std_items   <- as.integer(dif_items$standard     %||% integer(0))
      inter_items <- as.integer(dif_items$intersection %||% integer(0))
    } else {
      di      <- as.integer(dif_items)
      n_half  <- floor(length(di) / 2)
      std_items   <- di[seq_len(n_half)]
      inter_items <- di[seq(n_half + 1L, length(di))]
    }
  } else {
    # standard or intersection
    di <- as.integer(if (is.list(dif_items)) unlist(dif_items) else dif_items)
    std_items   <- if (dif_structure == "standard")     di else integer(0)
    inter_items <- if (dif_structure == "intersection") di else integer(0)
  }

  std_items   <- std_items[std_items   >= 1L & std_items   <= n_items]
  inter_items <- inter_items[inter_items >= 1L & inter_items <= n_items]
  all_dif     <- unique(c(std_items, inter_items))

  # --- Group assignment -------------------------------------------------------

  group     <- factor(sample(paste0("G", seq_len(n_groups)),
                              n_persons, replace = TRUE))
  group_int <- as.integer(group)
  is_focal  <- group_int > 1L

  # --- Demographic variables --------------------------------------------------
  # Created before response generation so the target cell is identifiable.

  demo_data <- data.frame(group = group)
  if (!is.null(demo_vars)) {
    for (v_name in names(demo_vars)) {
      demo_data[[v_name]] <- sample(demo_vars[[v_name]], n_persons, replace = TRUE)
    }
  }

  # --- Identify intersection target cell --------------------------------------

  if (dif_structure %in% c("intersection", "mixed") && !is.null(dif_group)) {
    in_target <- rep(TRUE, n_persons)
    for (v_name in names(dif_group)) {
      in_target <- in_target &
        as.character(demo_data[[v_name]]) == as.character(dif_group[[v_name]])
    }
  } else {
    in_target <- rep(FALSE, n_persons)
  }

  # --- Item parameters --------------------------------------------------------

  theta <- stats::rnorm(n_persons, mean = 0, sd = 1)
  a     <- stats::runif(n_items, 0.5, 2.0)
  b     <- stats::rnorm(n_items, 0, 1)

  a_mat <- matrix(a, nrow = n_persons, ncol = n_items, byrow = TRUE)
  b_mat <- matrix(b, nrow = n_persons, ncol = n_items, byrow = TRUE)

  # --- Apply standard DIF (focal groups) -------------------------------------

  for (j in std_items) {
    if (dif_type == "uniform") {
      b_mat[is_focal, j] <- b[j] + dif_effect
    } else {
      # Pure discrimination shift: multiply focal group's a-parameter only.
      # b is left unchanged so the ICCs cross, which is the defining
      # characteristic of non-uniform DIF.
      a_mat[is_focal, j] <- a[j] * (1 + dif_effect)
    }
  }

  # --- Apply intersection DIF (target cell only) -----------------------------

  for (j in inter_items) {
    if (dif_type == "uniform") {
      b_mat[in_target, j] <- b[j] + dif_effect
    } else {
      a_mat[in_target, j] <- a[j] * (1 + dif_effect)
    }
  }

  # --- Simulate responses ----------------------------------------------------

  theta_mat <- matrix(theta, nrow = n_persons, ncol = n_items)
  P         <- 1 / (1 + exp(-a_mat * (theta_mat - b_mat)))
  responses <- matrix(
    stats::rbinom(n_persons * n_items, 1, as.vector(P)),
    nrow = n_persons, ncol = n_items
  )

  # --- Build output data frame -----------------------------------------------

  df        <- as.data.frame(responses)
  names(df) <- paste0("item_", seq_len(n_items))
  df$group  <- group

  if (!is.null(demo_vars)) {
    for (v_name in names(demo_vars)) df[[v_name]] <- demo_data[[v_name]]
  }

  attr(df, "true_a")        <- a
  attr(df, "true_b")        <- b
  attr(df, "dif_items")     <- all_dif
  attr(df, "dif_effect")    <- dif_effect
  attr(df, "dif_type")      <- dif_type
  attr(df, "dif_structure") <- dif_structure
  attr(df, "dif_group")     <- dif_group
  attr(df, "demo_vars")     <- demo_vars

  df
}

# Null-coalescing helper (R < 4.4 doesn't have %||%)
`%||%` <- function(a, b) if (!is.null(a)) a else b

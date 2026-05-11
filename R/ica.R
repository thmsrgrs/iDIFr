#' Intersectional Contrast Analysis (ICA)
#'
#' @description
#' Runs a series of `idifr()` analyses — one per demographic variable and one
#' with the full intersectional formula — then classifies each item by
#' comparing where it was flagged.
#'
#' ICA distinguishes four patterns:
#' \describe{
#'   \item{`"amplified"`}{Flagged in both a single-variable run *and* the
#'     intersectional run. The intersectional design amplifies the signal.}
#'   \item{`"pure_intersection"`}{Flagged only in the intersectional run. DIF
#'     is invisible when examining demographic variables individually.}
#'   \item{`"obscured"`}{Flagged in a single-variable run but *not* in the
#'     intersectional run. The intersectional crossing dilutes or reverses the
#'     effect.}
#'   \item{`"none"`}{Not flagged in any analysis.}
#' }
#'
#' @param data A data frame containing item responses and demographic variables.
#' @param items Numeric or character vector identifying the item columns.
#' @param group A one-sided formula with **at least two** demographic variables,
#'   e.g. `~ gender * nationality * age_band`. ICA is not meaningful for a
#'   single-variable design — use [idifr()] directly in that case.
#' @param method `"LR"` or `"LRT"` — the DIF detection method passed to each
#'   `idifr()` call.
#' @param min_cell_size Minimum acceptable group size passed to `idifr()`.
#'   Default 50.
#' @param alpha Significance threshold. Default 0.05.
#' @param p_adjust P-value adjustment method passed to `idifr()`. Default
#'   `"BH"`.
#' @param effect_threshold Minimum effect size required for an item to be
#'   counted as flagged in a single-variable or intersectional run. For `"LR"`
#'   this is the minimum Nagelkerke ΔR²; for `"LRT"` the minimum std. chi.
#'   Defaults to `0.035` (ETS Class B). A stricter value here partially
#'   compensates for running multiple analyses without cross-analysis
#'   correction.
#' @param verbose Logical. Print progress and notes. Default `TRUE`.
#' @param ... Additional arguments passed to each `idifr()` call.
#'
#' @return An object of class `"idifr_ica"` containing:
#' \describe{
#'   \item{`classification`}{Data frame (one row per item) with columns
#'     `item`, `classification`, `marginal_variable`, `intersectional_flagged`,
#'     `method`.}
#'   \item{`single_results`}{Named list of `tidy()` data frames, one per
#'     demographic variable.}
#'   \item{`intersectional_result`}{`tidy()` data frame from the full
#'     intersectional run.}
#'   \item{`vars`}{Character vector of demographic variable names.}
#'   \item{`single_flags`}{Named list of character vectors: the item names
#'     flagged in each single-variable run (using `effect_threshold`).}
#'   \item{`method`, `alpha`, `effect_threshold`}{Parameters used.}
#'   \item{`n_groups`}{Number of intersectional groups.}
#'   \item{`call`}{The matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' library(iDIFr)
#' set.seed(42)
#' dat <- simulate_dif(1000, 20, 2, c(3, 7), 1.0, seed = 42)
#' dat$nationality <- sample(c("UK", "DE", "FR"), 1000, replace = TRUE)
#' dat$age_band    <- sample(c("Young", "Old"),   1000, replace = TRUE)
#'
#' ica_res <- ica(
#'   data   = dat,
#'   items  = 1:20,
#'   group  = ~ group * nationality * age_band,
#'   method = "LR"
#' )
#' print(ica_res)
#' tidy(ica_res)
#' }
#'
#' @seealso [idifr()] for single-run DIF analysis.
#' @export
ica <- function(data,
                items,
                group,
                method,
                min_cell_size    = 50,
                alpha            = 0.05,
                p_adjust         = "BH",
                effect_threshold = 0.035,
                verbose          = TRUE,
                ...) {

  call   <- match.call()
  method <- match.arg(method, c("LR", "LRT"))
  vars   <- all.vars(group)

  if (length(vars) < 2) {
    stop(
      "ICA requires at least 2 grouping variables.\n",
      "  The formula `~ ", deparse(group[[2]]), "` has only 1 variable.\n",
      "  Use idifr() directly for single-variable DIF analysis.",
      call. = FALSE
    )
  }

  if (verbose) {
    cli::cli_h1(paste0("ICA — Intersectional Contrast Analysis (", method, ")"))
    cli::cli_alert_info(
      "Running {length(vars) + 1L} idifr() analyses \\
       ({length(vars)} single-variable + 1 intersectional)."
    )
    cli::cli_alert_warning(
      "ICA runs multiple analyses without cross-analysis p-value correction. \\
       The {.arg effect_threshold} ({effect_threshold}) provides a de facto \\
       stricter criterion. Interpret 'pure_intersection' and 'obscured' \\
       classifications with caution in small samples."
    )
    cat("\n")
  }

  # --- Single-variable runs --------------------------------------------------

  single_runs <- lapply(vars, function(v) {
    f <- stats::as.formula(paste("~", v))
    if (verbose) cli::cli_alert_info("Single-variable run: {.field {v}}")
    idifr(data, items, f, method,
          min_cell_size = min_cell_size,
          alpha         = alpha,
          p_adjust      = p_adjust,
          verbose       = FALSE,
          ...)
  })
  names(single_runs) <- vars

  # --- Intersectional run ----------------------------------------------------

  if (verbose) cli::cli_alert_info("Running intersectional analysis...")
  inter_run <- idifr(data, items, group, method,
                     min_cell_size = min_cell_size,
                     alpha         = alpha,
                     p_adjust      = p_adjust,
                     verbose       = FALSE,
                     ...)
  if (verbose) cat("\n")

  # --- Extract tidy results --------------------------------------------------

  single_tidy <- lapply(single_runs, tidy)
  inter_tidy  <- tidy(inter_run)

  # Items from the intersectional run are the reference set
  all_items <- inter_run$items

  # --- Determine flagged items per run (honouring effect_threshold) ----------

  single_flags <- lapply(vars, function(v) {
    .ica_flagged_items(single_tidy[[v]], method, alpha, effect_threshold)
  })
  names(single_flags) <- vars

  inter_flagged <- .ica_flagged_items(inter_tidy, method, alpha, effect_threshold)

  # --- Classify items --------------------------------------------------------

  cls_rows <- lapply(all_items, function(it) {

    marginal      <- any(sapply(vars, function(v) it %in% single_flags[[v]]))
    intersectional <- it %in% inter_flagged

    classification <- dplyr::case_when(
      marginal  & intersectional  ~ "amplified",
      !marginal & intersectional  ~ "pure_intersection",
      marginal  & !intersectional ~ "obscured",
      TRUE                        ~ "none"
    )

    marginal_vars <- if (marginal) {
      paste(vars[sapply(vars, function(v) it %in% single_flags[[v]])],
            collapse = ", ")
    } else NA_character_

    data.frame(
      item                  = it,
      classification        = classification,
      marginal_variable     = marginal_vars,
      intersectional_flagged = intersectional,
      method                = method,
      stringsAsFactors = FALSE
    )
  })

  classification_df          <- do.call(rbind, cls_rows)
  rownames(classification_df) <- NULL

  structure(
    list(
      classification       = classification_df,
      single_results       = single_tidy,
      intersectional_result = inter_tidy,
      single_flags         = single_flags,
      vars                 = vars,
      method               = method,
      alpha                = alpha,
      effect_threshold     = effect_threshold,
      n_groups             = inter_run$groups$n_groups,
      call                 = call
    ),
    class = "idifr_ica"
  )
}


# Internal helper: item names that are flagged in a tidy result, using the
# caller-supplied effect_threshold (rather than the stored flagged column,
# which may use a different threshold).
.ica_flagged_items <- function(tidy_res, method, alpha, effect_threshold) {

  es_col <- switch(method,
    "LR"  = "delta_r2",
    "LRT" = "std_chi",
    NULL
  )

  if (is.null(es_col) || !es_col %in% names(tidy_res)) {
    return(character(0))
  }

  keep <- !is.na(tidy_res$p_adj)       & tidy_res$p_adj < alpha &
          !is.na(tidy_res[[es_col]])    & tidy_res[[es_col]] >= effect_threshold

  as.character(tidy_res$item[keep])
}


#' Print method for idifr_ica objects
#'
#' @param x An `idifr_ica` object.
#' @param ... Ignored.
#' @export
print.idifr_ica <- function(x, ...) {

  cli::cli_h1(paste0("ICA Results — ", x$method, " method"))

  cat(sprintf(
    "Design: %s  (%d intersectional groups)\n",
    paste(x$vars, collapse = " x "),
    x$n_groups
  ))
  cat(sprintf(
    "Alpha: %.2f   Effect threshold: %.3f   Adjustment: %s\n\n",
    x$alpha, x$effect_threshold,
    if (!is.null(x$call$p_adjust)) deparse(x$call$p_adjust) else "BH"
  ))

  cls_names  <- c("amplified", "pure_intersection", "obscured", "none")
  cls_labels <- c("Amplified", "Pure intersection", "Obscured", "None")

  w_cls  <- 20
  w_n    <- 8
  sep    <- strrep("─", w_cls + w_n + 40)

  cat(sprintf("%-*s  %-*s  %s\n", w_cls, "Classification", w_n, "Items",
              "Example items"))
  cat(sep, "\n")

  for (i in seq_along(cls_names)) {
    items_cls <- x$classification$item[
      x$classification$classification == cls_names[i]]
    n       <- length(items_cls)
    example <- if (n == 0) "—" else {
      ex <- paste(head(items_cls, 3), collapse = ", ")
      if (n > 3) paste0(ex, " …") else ex
    }
    cat(sprintf("%-*s  %-*d  %s\n", w_cls, cls_labels[i], w_n, n, example))
  }
  cat(sep, "\n")

  # Marginal DIF drivers
  marginal_rows <- x$classification[
    x$classification$classification %in% c("amplified", "obscured"), ]

  if (nrow(marginal_rows) > 0) {
    cat("\nMarginal DIF drivers:\n")
    for (v in x$vars) {
      v_items <- x$single_flags[[v]]
      # Restrict to items that are in the marginal rows
      reported <- intersect(v_items, marginal_rows$item)
      if (length(reported) > 0) {
        cat(sprintf("  %-16s %s\n",
                    paste0(v, ":"),
                    paste(reported, collapse = ", ")))
      }
    }
  }

  cat("\n")
  invisible(x)
}


#' Extract tidy classification table from an idifr_ica object
#'
#' @param x An `idifr_ica` object.
#' @param ... Ignored.
#' @return A data frame with one row per item.
#' @export
tidy.idifr_ica <- function(x, ...) {
  x$classification
}

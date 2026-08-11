# ============================================================================
# iDIFr Simulation Study — VERSION 4 (CLEAN)
#
# Single authoritative script covering all results sections:
#
#   sim_performance_v4.csv  → power / FPR for all conditions
#   sim_ica_v4.csv          → per-item ICA labels (intersectional conditions)
#   simulation_v4_tables.xlsx → all summary tables
#
# Conditions:
#   none         → Type I error (Section 4)
#   standard     → Main-effect power (Section 4.1)
#   intersection → Pure intersectional power (Section 4.2)
#   amplified    → Amplified power/accuracy (Section 4.3)
#   obscured     → Obscured power/accuracy (Section 4.4)
#
# Ground-truth ICA:
#   None              — null condition
#   Pure Intersectional — DIF confined to G1-UK cell only
#   Amplified         — main effect Δ for all G1 + cell boost Δ/2 for G1-UK
#   Obscured condition — even Δ for all G1; correct label = Amplified OR Obscured
#
# Seeds: condition_id * 10000L + rep_id  (fully reproducible)
# ============================================================================

library(iDIFr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(parallel)

# ── Design constants ───────────────────────────────────────────────────────
N_ITEMS   <- 40
DIF_ITEMS <- c(3, 7, 12, 18, 25, 30, 33, 37)
mag_levels   <- c("small", "medium", "large")
MAGNITUDES   <- c(small = 0.5, medium = 0.8, large = 1.2)
SAMPLE_SIZES <- c(500, 1000, 2000)
DIF_TYPES    <- c("uniform", "nonuniform")

# ============================================================================
# 1. HELPER: safe group size allocation (handles non-divisible n_persons)
# ============================================================================

safe_sizes <- function(n_persons, n_groups) {
  as.integer(diff(round(seq(0, n_persons, length.out = n_groups + 1))))
}

# ============================================================================
# 2. HELPER: generate 2PL responses
# ============================================================================

generate_responses <- function(theta, a_mat, b_mat) {
  n_p <- length(theta)
  n_i <- ncol(a_mat)
  P   <- 1 / (1 + exp(-a_mat * (matrix(theta, n_p, n_i) - b_mat)))
  matrix(rbinom(n_p * n_i, 1, as.vector(P)), nrow = n_p, ncol = n_i)
}

# ============================================================================
# 3. DATA GENERATING FUNCTIONS
# ============================================================================

# ── Null ────────────────────────────────────────────────────────────────────
generate_null <- function(n_persons, seed) {
  if (!is.null(seed)) set.seed(seed)
  sz          <- safe_sizes(n_persons, 4)
  group_ids   <- rep(1:4, times = sz)
  group       <- ifelse(group_ids %in% c(1, 2), "G1", "G2")
  nationality <- ifelse(group_ids %in% c(1, 3), "UK", "DE")
  theta <- rnorm(n_persons)
  a     <- runif(N_ITEMS, 0.5, 2.0)
  b     <- rnorm(N_ITEMS)
  a_mat <- matrix(a, n_persons, N_ITEMS, byrow = TRUE)
  b_mat <- matrix(b, n_persons, N_ITEMS, byrow = TRUE)
  resp  <- generate_responses(theta, a_mat, b_mat)
  df    <- as.data.frame(resp)
  names(df) <- paste0("item_", seq_len(N_ITEMS))
  df$group       <- group
  df$nationality <- nationality
  list(data = df, true_items = integer(0))
}

# ── Standard main-effect DIF — all group structures ─────────────────────────
generate_standard <- function(n_persons, magnitude, dif_type, seed,
                              group_structure = "intersectional") {
  if (!is.null(seed)) set.seed(seed)
  
  if (group_structure == "intersectional") {
    sz          <- safe_sizes(n_persons, 4)
    group_ids   <- rep(1:4, times = sz)
    group       <- ifelse(group_ids %in% c(1, 2), "G1", "G2")
    nationality <- ifelse(group_ids %in% c(1, 3), "UK", "DE")
  } else if (group_structure == "multigroup") {
    sz        <- safe_sizes(n_persons, 3)
    group_ids <- rep(1:3, times = sz)
    group     <- paste0("G", group_ids)
    nationality <- NULL
  } else {   # single
    sz        <- safe_sizes(n_persons, 2)
    group_ids <- rep(1:2, times = sz)
    group     <- paste0("G", group_ids)
    nationality <- NULL
  }
  
  theta <- rnorm(n_persons)
  a     <- runif(N_ITEMS, 0.5, 2.0)
  b     <- rnorm(N_ITEMS)
  a_mat <- matrix(a, n_persons, N_ITEMS, byrow = TRUE)
  b_mat <- matrix(b, n_persons, N_ITEMS, byrow = TRUE)
  
  is_G1 <- group == "G1"   # G1 is the focal group (harder items)
  
  for (j in DIF_ITEMS) {
    if (dif_type == "uniform") {
      b_mat[is_G1, j] <- b[j] + magnitude
    } else {
      a_mat[is_G1, j] <- a[j] * (1 + magnitude)
    }
  }
  
  resp      <- generate_responses(theta, a_mat, b_mat)
  df        <- as.data.frame(resp)
  names(df) <- paste0("item_", seq_len(N_ITEMS))
  df$group  <- group
  if (!is.null(nationality)) df$nationality <- nationality
  list(data = df, true_items = DIF_ITEMS)
}

# ── Pure intersectional — DIF confined to G1-UK only ────────────────────────
generate_intersection <- function(n_persons, magnitude, dif_type, seed) {
  if (!is.null(seed)) set.seed(seed)
  sz          <- safe_sizes(n_persons, 4)
  group_ids   <- rep(1:4, times = sz)
  group       <- ifelse(group_ids %in% c(1, 2), "G1", "G2")
  nationality <- ifelse(group_ids %in% c(1, 3), "UK", "DE")
  theta <- rnorm(n_persons)
  a     <- runif(N_ITEMS, 0.5, 2.0)
  b     <- rnorm(N_ITEMS)
  a_mat <- matrix(a, n_persons, N_ITEMS, byrow = TRUE)
  b_mat <- matrix(b, n_persons, N_ITEMS, byrow = TRUE)
  is_target <- group == "G1" & nationality == "UK"
  for (j in DIF_ITEMS) {
    if (dif_type == "uniform") {
      b_mat[is_target, j] <- b[j] + magnitude
    } else {
      a_mat[is_target, j] <- a[j] * (1 + magnitude)
    }
  }
  resp        <- generate_responses(theta, a_mat, b_mat)
  df          <- as.data.frame(resp)
  names(df)   <- paste0("item_", seq_len(N_ITEMS))
  df$group       <- group
  df$nationality <- nationality
  list(data = df, true_items = DIF_ITEMS)
}

# ── Amplified — main effect Δ for all G1 + cell boost Δ/2 for G1-UK ─────────
# Both marginal and intersectional analyses detect DIF.
# Intersectional adds resolution: G1-UK harder than G1-DE → Amplified
generate_amplified <- function(n_persons, magnitude, dif_type, seed) {
  if (!is.null(seed)) set.seed(seed)
  cell_boost  <- magnitude / 2
  sz          <- safe_sizes(n_persons, 4)
  group_ids   <- rep(1:4, times = sz)
  group       <- ifelse(group_ids %in% c(1, 2), "G1", "G2")
  nationality <- ifelse(group_ids %in% c(1, 3), "UK", "DE")
  theta <- rnorm(n_persons)
  a     <- runif(N_ITEMS, 0.5, 2.0)
  b     <- rnorm(N_ITEMS)
  a_mat <- matrix(a, n_persons, N_ITEMS, byrow = TRUE)
  b_mat <- matrix(b, n_persons, N_ITEMS, byrow = TRUE)
  is_G1    <- group == "G1"
  is_G1_UK <- group == "G1" & nationality == "UK"
  for (j in DIF_ITEMS) {
    if (dif_type == "uniform") {
      b_mat[is_G1, j]    <- b[j] + magnitude              # main effect
      b_mat[is_G1_UK, j] <- b[j] + magnitude + cell_boost # + cell boost
    } else {
      a_mat[is_G1, j]    <- a[j] * (1 + magnitude)
      a_mat[is_G1_UK, j] <- a[j] * (1 + magnitude + cell_boost)
    }
  }
  resp        <- generate_responses(theta, a_mat, b_mat)
  df          <- as.data.frame(resp)
  names(df)   <- paste0("item_", seq_len(N_ITEMS))
  df$group       <- group
  df$nationality <- nationality
  list(data = df, true_items = DIF_ITEMS)
}

# ── Obscured — even main effect Δ for all G1 (no cell boost) ─────────────────
# Marginal detects reliably (full G1 vs G2 sample).
# Intersectional may or may not confirm (smaller cells).
# At low power → Obscured; at high power → Amplified. Both correct.
generate_obscured <- function(n_persons, magnitude, dif_type, seed) {
  if (!is.null(seed)) set.seed(seed)
  sz          <- safe_sizes(n_persons, 4)
  group_ids   <- rep(1:4, times = sz)
  group       <- ifelse(group_ids %in% c(1, 2), "G1", "G2")
  nationality <- ifelse(group_ids %in% c(1, 3), "UK", "DE")
  theta <- rnorm(n_persons)
  a     <- runif(N_ITEMS, 0.5, 2.0)
  b     <- rnorm(N_ITEMS)
  a_mat <- matrix(a, n_persons, N_ITEMS, byrow = TRUE)
  b_mat <- matrix(b, n_persons, N_ITEMS, byrow = TRUE)
  is_G1 <- group == "G1"
  for (j in DIF_ITEMS) {
    if (dif_type == "uniform") {
      b_mat[is_G1, j] <- b[j] + magnitude   # even shift, no cell boost
    } else {
      a_mat[is_G1, j] <- a[j] * (1 + magnitude)
    }
  }
  resp        <- generate_responses(theta, a_mat, b_mat)
  df          <- as.data.frame(resp)
  names(df)   <- paste0("item_", seq_len(N_ITEMS))
  df$group       <- group
  df$nationality <- nationality
  list(data = df, true_items = DIF_ITEMS)
}

# ============================================================================
# 4. SINGLE REPLICATION
# ============================================================================

run_replication <- function(condition, rep_id) {
  
  seed_val <- condition$condition_id * 10000L + rep_id
  src      <- condition$dif_source
  n        <- condition$sample_size
  mag      <- condition$dif_magnitude
  dtype    <- condition$dif_type
  gstruct  <- condition$group_structure
  
  sim <- switch(src,
                none         = generate_null(n, seed_val),
                standard     = generate_standard(n, mag, dtype, seed_val, gstruct),
                intersection = generate_intersection(n, mag, dtype, seed_val),
                amplified    = generate_amplified(n, mag, dtype, seed_val),
                obscured     = generate_obscured(n, mag, dtype, seed_val),
                stop("Unknown dif_source: ", src)
  )
  
  dat       <- sim$data
  true_cols <- paste0("item_", sim$true_items)
  has_nat   <- "nationality" %in% names(dat)
  
  # Standard conditions: marginal formula, no ICA
  # All other intersectional conditions: full formula + ICA
  if (!has_nat || src == "standard") {
    grp_form <- ~ group
    run_ica  <- FALSE
  } else {
    grp_form <- ~ group * nationality
    run_ica  <- TRUE
  }
  
  result <- tryCatch(
    idifr(dat, items = 1:N_ITEMS, group = grp_form,
          method = c("LRT"), ica = run_ica, verbose = FALSE),
    error = function(e) e
  )
  
  if (inherits(result, "error")) {
    perf <- tibble(method = c("LR","LRT","MOB"),
                   power = NA, fpr = NA, concordance = NA,
                   error = conditionMessage(result))
    return(list(performance = perf, ica = NULL))
  }
  
  # ── Performance ──────────────────────────────────────────────────────────
  res          <- tidy(result)
  res$true_dif <- res$item %in% true_cols
  
  perf <- res %>%
    group_by(method) %>%
    summarise(
      power = if (sum(true_dif)  > 0) mean(flagged[true_dif])  else NA_real_,
      fpr   = if (sum(!true_dif) > 0) mean(flagged[!true_dif]) else NA_real_,
      .groups = "drop"
    )
  
  wide        <- res %>% select(item, method, flagged) %>%
    pivot_wider(names_from = method, values_from = flagged)
  method_cols <- intersect(c("LR","LRT","MOB"), names(wide))
  perf$concordance <- mean(apply(wide[method_cols], 1,
                                 function(x) length(unique(x)) == 1))
  perf$error <- NA_character_
  
  
  stats_out <- res %>%
    select(item, method, chi_sq, df, std_chi, chi_uniform, chi_nonuniform)
  
  
  # ── ICA output ───────────────────────────────────────────────────────────
  ica_out <- NULL
  if (run_ica) {
    ica_tbl <- tidy(result, table = "ica") %>%
      left_join(res %>% select(item, method, flagged, true_dif),
                by = c("item","method")) %>%
      mutate(
        item_idx     = as.integer(gsub("item_", "", item)),
        expected_ica = case_when(
          src == "none"         ~ "none",
          src == "intersection" ~ if_else(item_idx %in% DIF_ITEMS,
                                          "pure_intersection", "none"),
          src == "amplified"    ~ if_else(item_idx %in% DIF_ITEMS,
                                          "amplified", "none"),
          src == "obscured"     ~ if_else(item_idx %in% DIF_ITEMS,
                                          "amplified_or_obscured", "none"),
          TRUE ~ NA_character_
        ),
        correct = case_when(
          expected_ica == "amplified_or_obscured" ~
            ica_class %in% c("amplified","obscured"),
          TRUE ~ ica_class == expected_ica
        )
      )
    
    ica_out <- ica_tbl %>%
      select(item, item_idx, method, ica_class, expected_ica,
             true_dif, flagged, correct)
  }
  
  list(performance = perf, ica = ica_out, stats = stats_out)
}

# ============================================================================
# 5. CONDITION GRIDS
# ============================================================================

build_grids <- function() {
  
  null_grid <- expand.grid(
    group_structure  = c("single","multigroup","intersectional"),
    dif_source       = "none",
    magnitude_label  = "none",
    dif_type         = "uniform",
    sample_size      = SAMPLE_SIZES,
    stringsAsFactors = FALSE
  )
  null_grid$dif_magnitude <- NA_real_
  
  standard_grid <- expand.grid(
    group_structure  = c("single","multigroup","intersectional"),
    dif_source       = "standard",
    magnitude_label  = mag_levels,
    dif_type         = DIF_TYPES,
    sample_size      = SAMPLE_SIZES,
    stringsAsFactors = FALSE
  )
  standard_grid$dif_magnitude <- MAGNITUDES[standard_grid$magnitude_label]
  
  intersection_grid <- expand.grid(
    group_structure  = "intersectional",
    dif_source       = "intersection",
    magnitude_label  = mag_levels,
    dif_type         = DIF_TYPES,
    sample_size      = SAMPLE_SIZES,
    stringsAsFactors = FALSE
  )
  intersection_grid$dif_magnitude <- MAGNITUDES[intersection_grid$magnitude_label]
  
  amplified_grid <- expand.grid(
    group_structure  = "intersectional",
    dif_source       = "amplified",
    magnitude_label  = mag_levels,
    dif_type         = DIF_TYPES,
    sample_size      = SAMPLE_SIZES,
    stringsAsFactors = FALSE
  )
  amplified_grid$dif_magnitude <- MAGNITUDES[amplified_grid$magnitude_label]
  
  obscured_grid <- expand.grid(
    group_structure  = "intersectional",
    dif_source       = "obscured",
    magnitude_label  = mag_levels,
    dif_type         = DIF_TYPES,
    sample_size      = SAMPLE_SIZES,
    stringsAsFactors = FALSE
  )
  obscured_grid$dif_magnitude <- MAGNITUDES[obscured_grid$magnitude_label]
  
  grid <- bind_rows(null_grid, standard_grid, intersection_grid,
                    amplified_grid, obscured_grid)
  grid$condition_id <- seq_len(nrow(grid))
  grid
}

# ============================================================================
# 6. PARALLEL RUNNER
# ============================================================================

run_simulation_parallel <- function(grid, n_reps,
                                    perf_file, ica_file,
                                    n_cores    = NULL,
                                    batch_size = 200,
                                    resume     = FALSE) {
  
  if (is.null(n_cores)) n_cores <- max(1L, detectCores() - 1L)
  
  tasks <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    data.frame(condition_idx = i, rep = seq_len(n_reps),
               condition_id  = grid$condition_id[i])))
  
  done_keys <- character(0)
  if (resume && file.exists(perf_file)) {
    ex        <- read.csv(perf_file, stringsAsFactors = FALSE)
    done_keys <- unique(paste(ex$condition_id, ex$rep))
  }
  tasks$key <- paste(tasks$condition_id, tasks$rep)
  tasks     <- tasks[!(tasks$key %in% done_keys), ]
  
  if (nrow(tasks) == 0) {
    cat("All tasks already completed.\n")
    return(invisible(TRUE))
  }
  
  cat(sprintf("Running %d tasks on %d cores\n", nrow(tasks), n_cores))
  
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  
  clusterEvalQ(cl, { library(iDIFr); library(dplyr); library(tidyr) })
  clusterExport(cl, c("grid","run_replication",
                      "generate_null","generate_standard",
                      "generate_intersection","generate_amplified",
                      "generate_obscured","generate_responses",
                      "safe_sizes","N_ITEMS","DIF_ITEMS","MAGNITUDES"))
  
  meta_cols <- c("group_structure","dif_source","magnitude_label",
                 "dif_type","sample_size","dif_magnitude")
  n_batches <- ceiling(nrow(tasks) / batch_size)
  
  for (b in seq_len(n_batches)) {
    idx   <- ((b-1)*batch_size + 1):min(b*batch_size, nrow(tasks))
    batch <- tasks[idx, ]
    
    results <- parLapply(cl, seq_len(nrow(batch)), function(j) {
      cond <- grid[batch$condition_idx[j], ]
      run_replication(cond, batch$rep[j])
    })
    
    perf_rows <- do.call(rbind, lapply(seq_along(results), function(j) {
      cond <- grid[batch$condition_idx[j], ]
      bind_cols(cond[rep(1L, nrow(results[[j]]$performance)), meta_cols],
                results[[j]]$performance,
                condition_id = cond$condition_id, rep = batch$rep[j])
    }))
    write.table(perf_rows, perf_file, sep=",", row.names=FALSE,
                col.names=!file.exists(perf_file), append=file.exists(perf_file))
    
    ica_list <- Filter(Negate(is.null), lapply(seq_along(results), function(j) {
      if (is.null(results[[j]]$ica)) return(NULL)
      cond <- grid[batch$condition_idx[j], ]
      bind_cols(cond[rep(1L, nrow(results[[j]]$ica)), meta_cols],
                results[[j]]$ica,
                condition_id = cond$condition_id, rep = batch$rep[j])
    }))
    if (length(ica_list) > 0)
      write.table(do.call(rbind, ica_list), ica_file, sep=",",
                  row.names=FALSE, col.names=!file.exists(ica_file),
                  append=file.exists(ica_file))
    
    cat(sprintf("Batch %d / %d complete\n", b, n_batches))
  }
  invisible(TRUE)
}

# ============================================================================
# 7. RUN
# ============================================================================

grid <- build_grids()
cat(sprintf("Total conditions: %d\n", nrow(grid)))
cat("Conditions by source:\n"); print(table(grid$dif_source))

for (f in c("sim_performance_v4.csv","sim_ica_v4.csv")) {
  if (file.exists(f)) file.remove(f)
}

run_simulation_parallel(
  grid,
  n_reps     = 30,
  perf_file  = "sim_performance_v4_LRT_only.csv",
  ica_file   = "sim_ica_v4_LRT_only.csv",
  n_cores    = NULL,
  batch_size = 200,
  resume     = FALSE
)

# ============================================================================
# 8. AGGREGATION
# ============================================================================

perf <- read.csv("sim_performance_v4.csv", stringsAsFactors = FALSE)
ica  <- read.csv("sim_ica_v4.csv",         stringsAsFactors = FALSE)

mf <- function(x) factor(x, levels = mag_levels)

recode_ica <- function(x) recode(x,
                                 "none"              = "None",
                                 "amplified"         = "Amplified",
                                 "pure_intersection" = "Pure Intersectional",
                                 "obscured"          = "Obscured")

ica_dist <- function(data, src) {
  data %>%
    filter(dif_source == src, true_dif == TRUE) %>%
    mutate(magnitude_label = mf(magnitude_label),
           ica_class = recode_ica(ica_class)) %>%
    group_by(dif_type, magnitude_label, sample_size, method, ica_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(dif_type, magnitude_label, sample_size, method) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    select(-n) %>%
    pivot_wider(names_from = ica_class, values_from = pct, values_fill = 0) %>%
    arrange(dif_type, magnitude_label, sample_size, method)
}

ica_dist_summary <- function(data, src) {
  data %>%
    filter(dif_source == src, true_dif == TRUE) %>%
    mutate(magnitude_label = mf(magnitude_label),
           ica_class = recode_ica(ica_class)) %>%
    group_by(dif_type, magnitude_label, method, ica_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(dif_type, magnitude_label, method) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    select(-n) %>%
    pivot_wider(names_from = ica_class, values_from = pct, values_fill = 0) %>%
    arrange(dif_type, magnitude_label, method)
}

# Type I error
type1 <- perf %>%
  filter(dif_source == "none") %>%
  group_by(group_structure, sample_size, method) %>%
  summarise(mean_fpr = round(mean(fpr, na.rm=TRUE), 4),
            sd_fpr   = round(sd(fpr,   na.rm=TRUE), 4),
            .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = c(mean_fpr, sd_fpr))

# Section 4.1: main-effect power
standard_power <- perf %>%
  filter(dif_source == "standard") %>%
  mutate(magnitude_label = mf(magnitude_label)) %>%
  group_by(group_structure, dif_type, magnitude_label, method) %>%
  summarise(mean_power = round(mean(power, na.rm=TRUE), 3), .groups="drop") %>%
  pivot_wider(names_from = method, values_from = mean_power) %>%
  arrange(group_structure, dif_type, magnitude_label)

# Section 4.2: pure intersectional power
intersection_power <- perf %>%
  filter(dif_source == "intersection") %>%
  mutate(magnitude_label = mf(magnitude_label)) %>%
  group_by(dif_type, magnitude_label, sample_size, method) %>%
  summarise(mean_power = round(mean(power, na.rm=TRUE), 3), .groups="drop") %>%
  pivot_wider(names_from = method, values_from = mean_power) %>%
  arrange(dif_type, magnitude_label, sample_size)

# Section 4.3: amplified power + ICA
amplified_power <- perf %>%
  filter(dif_source == "amplified") %>%
  mutate(magnitude_label = mf(magnitude_label)) %>%
  group_by(dif_type, magnitude_label, sample_size, method) %>%
  summarise(mean_power = round(mean(power, na.rm=TRUE), 3), .groups="drop") %>%
  pivot_wider(names_from = method, values_from = mean_power) %>%
  arrange(dif_type, magnitude_label, sample_size)

amplified_ica         <- ica_dist(ica, "amplified")
amplified_ica_summary <- ica_dist_summary(ica, "amplified")

# Section 4.4: obscured power + ICA
obscured_power <- perf %>%
  filter(dif_source == "obscured") %>%
  mutate(magnitude_label = mf(magnitude_label)) %>%
  group_by(dif_type, magnitude_label, sample_size, method) %>%
  summarise(mean_power = round(mean(power, na.rm=TRUE), 3), .groups="drop") %>%
  pivot_wider(names_from = method, values_from = mean_power) %>%
  arrange(dif_type, magnitude_label, sample_size)

obscured_ica         <- ica_dist(ica, "obscured")
obscured_ica_summary <- ica_dist_summary(ica, "obscured")

# Print all
cat("\n=== TYPE I ERROR ===\n");                    print(type1,                  n=100)
cat("\n=== 4.1 MAIN-EFFECT POWER ===\n");          print(standard_power,         n=100)
cat("\n=== 4.2 PURE INTERSECTIONAL POWER ===\n");  print(intersection_power,     n=100)
cat("\n=== 4.3 AMPLIFIED POWER ===\n");            print(amplified_power,        n=100)
cat("\n=== 4.3 AMPLIFIED ICA (by N) ===\n");       print(amplified_ica,          n=200)
cat("\n=== 4.3 AMPLIFIED ICA (summary) ===\n");    print(amplified_ica_summary,  n=100)
cat("\n=== 4.4 OBSCURED POWER ===\n");             print(obscured_power,         n=100)
cat("\n=== 4.4 OBSCURED ICA (by N) ===\n");        print(obscured_ica,           n=200)
cat("\n=== 4.4 OBSCURED ICA (summary) ===\n");     print(obscured_ica_summary,   n=100)

# ============================================================================
# 9. EXPORT
# ============================================================================

wb <- createWorkbook()
add_tbl <- function(wb, sheet, df) {
  addWorksheet(wb, sheet)
  writeDataTable(wb, sheet, df)
  setColWidths(wb, sheet, cols=seq_along(df), widths="auto")
}

add_tbl(wb, "Type I Error",           type1)
add_tbl(wb, "4.1 Main Effect Power",  standard_power)
add_tbl(wb, "4.2 Intersection Power", intersection_power)
add_tbl(wb, "4.3 Amplified Power",    amplified_power)
add_tbl(wb, "4.3 Amplified ICA",      amplified_ica)
add_tbl(wb, "4.3 Amplified Summary",  amplified_ica_summary)
add_tbl(wb, "4.4 Obscured Power",     obscured_power)
add_tbl(wb, "4.4 Obscured ICA",       obscured_ica)
add_tbl(wb, "4.4 Obscured Summary",   obscured_ica_summary)

saveWorkbook(wb, "simulation_v4_tables.xlsx", overwrite = TRUE)
cat("\nDone. Files: sim_performance_v4.csv  sim_ica_v4.csv",
    " simulation_v4_tables.xlsx\n")
# validate_mob_nonuniform.R
# Validation suite for the MOB dual-test non-uniform DIF fix.
# Sprint: parallel raw + absolute CUSUM tests with DIF-type classification.
#
# Run with: source("validate_mob_nonuniform.R")

devtools::load_all(".")

cat("=================================================================\n")
cat("MOB dual-test fix -- 5 validation tests\n")
cat("=================================================================\n\n")


# ---- Shared datasets --------------------------------------------------------

# dat_std: 1000 persons, 20 items, UNIFORM DIF on items 3, 7, 12
dat_std <- simulate_dif(n_persons  = 1000,
                        n_items    = 20,
                        n_groups   = 2,
                        dif_items  = c(3, 7, 12),
                        dif_effect = 1.0,
                        dif_type   = "uniform",
                        seed       = 42)

# dat_nu: 1000 persons, 20 items, NON-UNIFORM DIF on items 3, 7
dat_nu <- simulate_dif(n_persons  = 1000,
                       n_items    = 20,
                       n_groups   = 2,
                       dif_items  = c(3, 7),
                       dif_effect = 1.0,
                       dif_type   = "nonuniform",
                       seed       = 42)

# dat_mixed: 50 items, MIXED DIF
#   Uniform:          items 5, 15, 25  (b-shift only)
#   Non-uniform:      items 10, 35     (a-shift only -- crossing ICC)
#   Both:             items 20, 40     (a-shift + b-shift)
#   Null:             all others
set.seed(99)
n     <- 1000
theta <- rnorm(n)
grp   <- sample(c("G1", "G2"), n, replace = TRUE)
a0    <- runif(50, 0.7, 1.5)
b0    <- rnorm(50, 0, 0.8)

# Person-level parameter matrices (group G1 = reference, G2 = focal)
is_focal <- grp == "G2"
a_mat <- matrix(a0, nrow = n, ncol = 50, byrow = TRUE)
b_mat <- matrix(b0, nrow = n, ncol = 50, byrow = TRUE)

# Uniform DIF: b-shift for focal group
for (j in c(5, 15, 25)) b_mat[is_focal, j] <- b0[j] + 1.2

# Non-uniform DIF: a-shift only (no b-shift → ICC crosses at b0)
for (j in c(10, 35)) a_mat[is_focal, j] <- a0[j] * 2.5

# Both: a-shift + b-shift
for (j in c(20, 40)) {
  a_mat[is_focal, j] <- a0[j] * 2.5
  b_mat[is_focal, j] <- b0[j] + 0.8
}

P          <- 1 / (1 + exp(-(a_mat * (matrix(theta, n, 50) - b_mat))))
resp       <- matrix(rbinom(n * 50, 1, as.vector(P)), nrow = n, ncol = 50)
colnames(resp) <- paste0("item_", 1:50)
dat_mixed  <- as.data.frame(resp)
dat_mixed$group <- grp

# dat_null: 1000 persons, 20 items, NO DIF
dat_null <- simulate_dif(n_persons  = 1000,
                         n_items    = 20,
                         n_groups   = 2,
                         dif_items  = integer(0),
                         dif_effect = 0,
                         seed       = 99)


# =============================================================================
# Test 1: Uniform DIF
# =============================================================================
cat("--- Test 1: Uniform DIF ---\n")

result_mob_u <- idifr(dat_std, paste0("item_", 1:20), ~ group,
                      "MOB", verbose = TRUE)
cat("\n")

mob_u <- tidy(result_mob_u)
flagged_u <- as.character(mob_u$item[!is.na(mob_u$flagged) & mob_u$flagged])
cat("Flagged:", paste(flagged_u, collapse = ", "), "\n")
cat("True:    item_3, item_7, item_12\n")

correct_u <- all(c("item_3", "item_7", "item_12") %in% flagged_u)
n_fp      <- length(setdiff(flagged_u, c("item_3", "item_7", "item_12")))
# At std_diff >= 0.20 (BH α = 0.05, 20 items) ≤1 borderline FP is expected
# and accepted; will be characterised in the simulation study.
fp_ok     <- n_fp <= 1L
cat("Uniform flagged correctly:  ", if (correct_u) "[PASS]" else "[FAIL]", "\n")
cat("False positives (≤1 OK):    ", if (fp_ok) "[PASS]" else "[FAIL]",
    sprintf("(%d FP)", n_fp), "\n")

# Check DIF type for flagged items is Uniform (or Uniform and Non-uniform)
types_u <- mob_u$dif_type[mob_u$item %in% c("item_3", "item_7", "item_12")]
all_uniform_type <- all(types_u %in% c("Uniform", "Uniform and Non-uniform"),
                        na.rm = TRUE)
cat("Uniform items have correct dif_type:", if (all_uniform_type) "[PASS]" else "[FAIL]",
    "(types:", paste(types_u, collapse = ", "), ")\n")
cat("\n")


# =============================================================================
# Test 2: Non-uniform DIF
# =============================================================================
cat("--- Test 2: Non-uniform DIF ---\n")

result_mob_nu <- idifr(dat_nu, paste0("item_", 1:20), ~ group,
                       "MOB", verbose = TRUE)
cat("\n")

mob_nu <- tidy(result_mob_nu)
flagged_nu <- as.character(mob_nu$item[!is.na(mob_nu$flagged) & mob_nu$flagged])
cat("Flagged:", paste(flagged_nu, collapse = ", "), "\n")
cat("True:    item_3, item_7\n")

correct_nu <- all(c("item_3", "item_7") %in% flagged_nu)
cat("Non-uniform flagged:        ", if (correct_nu) "[PASS]" else "[FAIL]", "\n")

# Check DIF type is Non-uniform
types_nu <- mob_nu$dif_type[mob_nu$item %in% c("item_3", "item_7")]
nu_type_ok <- all(types_nu %in% c("Non-uniform", "Uniform and Non-uniform"),
                  na.rm = TRUE)
cat("Non-uniform dif_type:       ", if (nu_type_ok) "[PASS]" else "[FAIL]",
    "(types:", paste(types_nu, collapse = ", "), ")\n")
cat("\n")


# =============================================================================
# Test 3: Mixed DIF
# =============================================================================
cat("--- Test 3: Mixed DIF (50 items) ---\n")

result_mob_mixed <- idifr(dat_mixed, paste0("item_", 1:50), ~ group,
                           "MOB", verbose = FALSE)
mob_mixed <- tidy(result_mob_mixed)
mob_flagged <- mob_mixed[!is.na(mob_mixed$flagged) & mob_mixed$flagged, ]

cat("\nFlagged items:\n")
print(mob_flagged[, c("item", "dif_type", "std_diff", "p_adj")])

cat("\nTrue Uniform:     item_5, item_15, item_25\n")
cat("True Non-uniform: item_10, item_35\n")
cat("True Both:        item_20, item_40\n\n")

flagged_mixed <- as.character(mob_flagged$item)

uni_ok  <- all(c("item_5", "item_15", "item_25") %in% flagged_mixed)
nu_ok   <- all(c("item_10", "item_35") %in% flagged_mixed)
both_ok <- all(c("item_20", "item_40") %in% flagged_mixed)

cat("Uniform items flagged:      ", if (uni_ok) "[PASS]" else "[FAIL]", "\n")
cat("Non-uniform items flagged:  ", if (nu_ok)  "[PASS]" else "[FAIL]", "\n")
cat("Both items flagged:         ", if (both_ok) "[PASS]" else "[FAIL]", "\n")

# Check type classifications
t5  <- mob_mixed$dif_type[mob_mixed$item == "item_5"]
t10 <- mob_mixed$dif_type[mob_mixed$item == "item_10"]
t20 <- mob_mixed$dif_type[mob_mixed$item == "item_20"]

cat("item_5  type (expect Uniform):")
cat(" '", t5, "'", if (t5 %in% c("Uniform","Uniform and Non-uniform")) "[PASS]" else "[FAIL]", "\n")
cat("item_10 type (expect Non-uniform):")
cat(" '", t10, "'", if (t10 %in% c("Non-uniform","Uniform and Non-uniform")) "[PASS]" else "[FAIL]", "\n")
cat("item_20 type (expect Uniform and Non-uniform):")
cat(" '", t20, "'", if (t20 == "Uniform and Non-uniform") "[PASS]" else "[FAIL]", "\n")
cat("\n")


# =============================================================================
# Test 4: Null -- no DIF
# =============================================================================
cat("--- Test 4: Null (no DIF) ---\n")

result_mob_null <- idifr(dat_null, paste0("item_", 1:20), ~ group,
                          "MOB", verbose = FALSE)
mob_null    <- tidy(result_mob_null)
flagged_null <- as.character(mob_null$item[!is.na(mob_null$flagged) & mob_null$flagged])

cat("Null flagged:", if (length(flagged_null) == 0) "none" else paste(flagged_null, collapse = ", "), "\n")
cat("Should be:    none\n")
cat("Type I error free:", if (length(flagged_null) == 0) "[PASS]" else "[FAIL]", "\n\n")


# =============================================================================
# Test 5: Direction tables
# =============================================================================
cat("--- Test 5: Direction tables ---\n\n")
cat("--- Uniform DIF result ---\n")
print(result_mob_u)

cat("\n--- Non-uniform DIF result ---\n")
print(result_mob_nu)

cat("\n=================================================================\n")
cat("Validation complete.\n")
cat("=================================================================\n")

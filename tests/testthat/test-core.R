test_that("simulate_dif produces correct structure", {
  dat <- simulate_dif(n_persons = 200, n_items = 10, n_groups = 2, seed = 42)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 200)
  expect_true("group" %in% names(dat))
  expect_true(all(paste0("item_", 1:10) %in% names(dat)))
  expect_true(all(unlist(dat[paste0("item_", 1:10)]) %in% c(0, 1, NA)))
})

test_that(".build_groups handles single variable correctly", {
  dat <- simulate_dif(100, 10, 2, seed = 1)
  grps <- iDIFr:::.build_groups(dat, ~ group, min_cell_size = 20)
  expect_s3_class(grps, "idifr_groups")
  expect_equal(grps$vars, "group")
  expect_equal(grps$n_groups, 2)
})

test_that(".build_groups handles intersectional formula", {
  dat <- simulate_dif(300, 10, 2, seed = 2)
  dat$nationality <- sample(c("UK", "DE", "FR"), 300, replace = TRUE)
  grps <- iDIFr:::.build_groups(dat, ~ group * nationality,
                                 min_cell_size = 20)
  expect_equal(grps$vars, c("group", "nationality"))
  expect_lte(grps$n_groups, 6)
})

test_that(".build_groups flags small cells", {
  dat <- simulate_dif(100, 5, 2, seed = 3)
  dat$rare <- c(rep("A", 95), rep("B", 5))
  grps <- iDIFr:::.build_groups(dat, ~ rare, min_cell_size = 50)
  expect_true(grps$has_small_cells)
})

test_that(".resolve_items accepts numeric indices", {
  dat <- simulate_dif(100, 10, seed = 1)
  cols <- iDIFr:::.resolve_items(dat, 1:5)
  expect_equal(cols, paste0("item_", 1:5))
})

test_that(".resolve_items accepts character names", {
  dat <- simulate_dif(100, 10, seed = 1)
  cols <- iDIFr:::.resolve_items(dat, c("item_1", "item_3"))
  expect_equal(cols, c("item_1", "item_3"))
})

test_that(".resolve_items rejects non-binary items", {
  dat <- simulate_dif(100, 5, seed = 1)
  dat$item_1 <- dat$item_1 + 2
  expect_error(iDIFr:::.resolve_items(dat, 1:5), "dichotomously scored")
})

test_that("idifr errors without method", {
  dat <- simulate_dif(200, 10, seed = 1)
  expect_error(
    idifr(dat, items = 1:10, group = ~ group),
    "must specify"
  )
})

test_that("idifr errors with unknown method", {
  dat <- simulate_dif(200, 10, seed = 1)
  expect_error(
    idifr(dat, items = 1:10, group = ~ group, method = "MH"),
    "Unknown method"
  )
})

test_that("idifr runs LR and returns idifr object", {
  dat <- simulate_dif(400, 10, dif_items = 3, dif_effect = 1.2, seed = 42)
  result <- idifr(dat, items = 1:10, group = ~ group,
                  method = "LR", verbose = FALSE)
  expect_s3_class(result, "idifr")
  expect_true("results" %in% names(result))
  expect_true("groups"  %in% names(result))
  expect_equal(nrow(result$results), 10)
})

test_that("idifr runs LRT and returns idifr object", {
  dat <- simulate_dif(400, 10, dif_items = 3, dif_effect = 1.2, seed = 42)
  result <- idifr(dat, items = 1:10, group = ~ group,
                  method = "LRT", verbose = FALSE)
  expect_s3_class(result, "idifr")
  expect_true("results" %in% names(result))
  expect_equal(nrow(result$results), 10)
  expect_true(all(c("chi_sq", "std_chi", "es_class", "flagged") %in%
                    names(result$results)))
})

test_that("tidy.idifr returns a data frame", {
  dat <- simulate_dif(400, 10, seed = 1)
  result <- idifr(dat, items = 1:10, group = ~ group,
                  method = "LR", verbose = FALSE)
  out <- tidy(result)
  expect_s3_class(out, "data.frame")
  expect_true("item" %in% names(out))
  expect_true("method" %in% names(out))
  expect_true("flagged" %in% names(out))
})

test_that("check_groups returns idifr_groups invisibly", {
  dat <- simulate_dif(300, 10, 2, seed = 5)
  grp <- check_groups(dat, group = ~ group, plot = FALSE)
  expect_s3_class(grp, "idifr_groups")
})

test_that("merge_groups recodes levels correctly", {
  dat <- simulate_dif(300, 10, 2, seed = 6)
  dat$age <- sample(c("18-24", "25-30", "31-45"), 300, replace = TRUE)
  grp <- check_groups(dat, ~ age, plot = FALSE)
  merged <- merge_groups(grp, age = list("18-30" = c("18-24", "25-30")))
  expect_true("18-30" %in% levels(merged$age))
  expect_false("18-24" %in% levels(merged$age))
})

test_that("MOB uniform direction labels align with residual sign", {
  tbl <- iDIFr:::.mob_direction_table(
    item_name  = "item_1",
    scores_raw = c(0.2, 0.1, -0.1, -0.2),
    var_data   = c("A", "A", "B", "B"),
    depth      = 1,
    dif_type   = "Uniform"
  )

  expect_equal(tbl$direction[tbl$group == "A"], "Advantaged")
  expect_equal(tbl$direction[tbl$group == "B"], "Disadvantaged")
})

toy_binary <- function() correct_metaweb_proba(toy_metaweb(), threshold = 0.5)

test_that("null model output is well formed and reproducible", {
  pa <- toy_sites()
  bin <- toy_binary()

  set.seed(1)
  r1 <- get_indic_null(pa, bin, n_null = 20)
  set.seed(1)
  r2 <- get_indic_null(pa, bin, n_null = 20)

  expect_identical(r1, r2)
  expect_named(r1, c("site", "index", "observed", "null_mean", "null_sd", "ses", "p_value"))
  expect_setequal(r1$site, c("s1", "s2", "s3"))
  expect_false(any(c("Species", "Species_taxa", "Link_max") %in% r1$index))
  expect_true(all(r1$p_value > 0 & r1$p_value <= 1, na.rm = TRUE))
})

test_that("null communities keep the richness of each site", {
  set.seed(2)
  for (method in c("equiprobable", "frequency")) {
    r <- get_indic_null(toy_sites(), toy_binary(), method = method, n_null = 10, indices = "Species")
    expect_equal(r$null_mean, r$observed)
    expect_true(all(is.na(r$ses)))
  }
})

test_that("get_indic_cells runs the null model only when asked", {
  pa <- toy_sites()
  bin <- toy_binary()
  expect_s3_class(get_indic_cells(pa, bin), "data.frame")

  set.seed(3)
  res <- get_indic_cells(pa, bin, null_model = "frequency", n_null = 10)
  expect_named(res, c("indices", "null_model"))
  expect_identical(res$indices, get_indic_cells(pa, bin))
  expect_s3_class(plot_indic_null(res$null_model, indices = c("Connectance", "TL_moy")), "ggplot")

  expect_error(get_indic_cells(pa, bin, null_model = "curveball"), "should be one of")
  expect_error(get_indic_null(pa, bin, indices = "not_an_index"), "Unknown index")
})

test_that("SES and p-values are computed correctly", {
  s <- summarise_null("site", c(x = 10), cbind(x = 1:9))
  expect_equal(s$null_mean, 5)
  expect_equal(s$ses, (10 - 5) / sd(1:9))
  expect_equal(s$p_value, 2 * (0 + 1) / (9 + 1))
})

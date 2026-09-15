test_that("apply_model_metaweb returns probabilities with prey in rows", {
  tr <- toy_traits()
  mw <- apply_model_metaweb(tr, toy_pars)

  expect_equal(dim(mw), c(8, 8))
  expect_identical(rownames(mw), rownames(tr))
  expect_identical(colnames(mw), rownames(tr))
  expect_true(all(mw >= 0 & mw <= 1))
  # cod (100 cm) is far more likely to eat capelin (15 cm) than the reverse
  expect_gt(mw["Mallotus_villosus", "Gadus_morhua"], mw["Gadus_morhua", "Mallotus_villosus"])

  o <- toy_pars[["a0"]] + toy_pars[["a1"]] * log10(100)
  r <- toy_pars[["b0"]] + toy_pars[["b1"]] * log10(100)
  expect_equal(mw["Mallotus_villosus", "Gadus_morhua"], exp(-(o - log10(15))^2 / (2 * r^2)))
})

test_that("parameters can be a vector, a file or the calibration output", {
  tr <- toy_traits()
  expected <- apply_model_metaweb(tr, toy_pars)

  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path))
  write.table(toy_pars, path)
  expect_equal(apply_model_metaweb(tr, path), expected)
  expect_equal(apply_model_metaweb(tr, list(calibrated_data = toy_pars)), expected)
})

test_that("trophic levels of a small web are correct", {
  nm <- c("P", "A", "B", "C")
  w <- matrix(0, 4, 4, dimnames = list(nm, nm))
  w["P", "A"] <- 1
  w["A", "B"] <- 1
  w["A", "C"] <- 1
  w["B", "C"] <- 1
  expect_equal(trophic_levels(w)$TL, c(1, 2, 3, 3.5))
})

test_that("trophic_levels() matches NetIndices::TrophInd()", {
  skip_if_not_installed("NetIndices")
  set.seed(42)
  for (k in 1:20) {
    n <- sample(5:40, 1)
    w <- (matrix(runif(n * n), n) > 0.8) * 1
    diag(w)[1:2] <- 1
    dimnames(w) <- list(paste0("s", 1:n), paste0("s", 1:n))
    ref <- NetIndices::TrophInd(Tij = t(w))
    res <- trophic_levels(w)
    expect_equal(res$TL, ref$TL, tolerance = 1e-8)
    expect_equal(res$OI, ref$OI, tolerance = 1e-8)
  }
})

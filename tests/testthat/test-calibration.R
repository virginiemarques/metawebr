test_that("calibration data are deduplicated and down-sampled", {
  d <- data.frame(
    predator = rep(c("P1", "P2"), c(12, 3)),
    prey = "x",
    standardised_predator_length = c(rep(50, 12), 60, 60, 60),
    si_prey_length = c(1:12, 5, 5, 6)
  )
  set.seed(1)
  out <- prepare_calibration_data(d, observation_resample = 5)

  # P1 capped at 5 records; P2 has one duplicate
  expect_equal(out$MPred, log10(c(rep(50, 5), 60, 60)))
  expect_error(prepare_calibration_data(d[, -1]), "predator")
})

test_that("model() equals the niche-model likelihood", {
  d <- data.frame(MPrey = log10(c(5, 10, 20, 8)), MPred = log10(c(50, 80, 200, 60)))
  p <- c(-0.5, 0.9, 0.2, 0.1)

  o <- p[1] + p[2] * d$MPred
  r <- p[3] + p[4] * d$MPred
  m <- mean(unique(d$MPrey))
  s <- sd(unique(d$MPrey))
  pML <- exp(-(o - d$MPrey)^2 / 2 / r^2) * dnorm(d$MPrey, m, s) /
    (r / (r^2 + s^2)^0.5 * exp(-(o - m)^2 / 2 / (r^2 + s^2)))

  expect_equal(model(p, d), -sum(log(pML)))
})

test_that("metaweb_mod_parameters output can be used by apply_model_metaweb", {
  skip_on_cran()
  set.seed(1)
  res <- metaweb_mod_parameters(max_time = 2, verbose = FALSE)

  expect_named(res$calibrated_data, c("a0", "a1", "b0", "b1"))
  expect_true(is.finite(res$neg_log_likelihood))
  expect_true(is.matrix(apply_model_metaweb(toy_traits(), res)))
})

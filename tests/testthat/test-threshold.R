test_that("get_proba returns a threshold from the grid", {
  thr <- get_proba(toy_metaweb(), toy_traits(), plot_diag = FALSE)
  expect_length(thr, 1)
  expect_true(thr %in% c(seq(0.1, 0.975, 0.025), 0.99, 0.995))
})

test_that("get_proba matches taxa by name (regression)", {
  tr <- toy_traits()
  mw <- toy_metaweb()
  expect_identical(get_proba(mw, tr, plot_diag = FALSE),
                   get_proba(mw, tr[rev(rownames(tr)), ], plot_diag = FALSE))
})

test_that("correct_metaweb_proba binarises at the chosen threshold", {
  mw <- toy_metaweb()
  thr <- get_proba(mw, toy_traits(), plot_diag = FALSE)

  # regression: the message used to round the threshold (0.675 shown as 0.7)
  expect_message(bin <- correct_metaweb_proba(mw, toy_traits(), plot_diag = FALSE),
                 paste("Threshold value chosen is", thr), fixed = TRUE)
  attr(mw, "corrections") <- NULL
  expect_identical(bin, (mw >= thr) * 1)

  expect_identical(correct_metaweb_proba(mw, threshold = 0.5), (mw >= 0.5) * 1)
  expect_error(correct_metaweb_proba(mw, threshold = 2), "between 0 and 1")
})

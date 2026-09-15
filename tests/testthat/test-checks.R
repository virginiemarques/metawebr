test_that("traits without taxon row names are rejected", {
  tr <- toy_traits()
  rownames(tr) <- NULL
  expect_error(apply_model_metaweb(tr, toy_pars), "row names")

  tr <- toy_traits()
  rownames(tr) <- as.character(seq_len(nrow(tr)))
  expect_error(apply_model_metaweb(tr, toy_pars), "row names")
  expect_error(get_proba(toy_metaweb(), tr), "row names")
})

test_that("missing columns and missing lengths are reported", {
  expect_error(apply_model_metaweb(toy_traits()[, c("species", "TrophicLevel")], toy_pars),
               "missing required column")

  tr <- toy_traits()
  tr["Clupea_harengus", "CommonLengthEstim"] <- NA
  expect_error(apply_model_metaweb(tr, toy_pars), "Clupea_harengus")
})

test_that("malformed metawebs are rejected", {
  mw <- apply_model_metaweb(toy_traits(), toy_pars)
  expect_error(correct_metaweb_fish(mw[, -1], toy_traits()), "square")

  shuffled <- mw
  colnames(shuffled) <- rev(colnames(shuffled))
  expect_error(get_proba(shuffled, toy_traits()), "identical")

  with_na <- mw
  with_na[1, 2] <- NA
  expect_error(plot_tree_network(with_na), "missing values")

  expect_error(correct_metaweb_proba(unname(mw), threshold = 0.5), "names")
})

test_that("model parameters are validated", {
  expect_error(apply_model_metaweb(toy_traits(), c(1, 2, 3)), "four numbers")
  expect_error(apply_model_metaweb(toy_traits(), "does/not/exist.txt"), "does not exist")
})

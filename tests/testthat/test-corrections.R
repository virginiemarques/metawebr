only_correction <- function(mw, tr, ...) {
  defaults <- list(herbivory = FALSE, water_pos = FALSE, small_fish = FALSE,
                   PP_SP_add = FALSE, exception_BigFish = FALSE)
  args <- utils::modifyList(defaults, list(...))
  suppressMessages(do.call(correct_metaweb_fish, c(list(mw, tr), args)))
}

test_that("herbivores lose their prey (regression: correction was never applied)", {
  tr <- toy_traits()
  tr["Boreogadus_saida", "TrophicLevel"] <- 2
  mw <- apply_model_metaweb(tr, toy_pars)
  expect_gt(sum(mw[, "Boreogadus_saida"]), 0)

  out <- only_correction(mw, tr, herbivory = TRUE)
  expect_equal(sum(out[, "Boreogadus_saida"]), 0)
  expect_equal(out[, "Gadus_morhua"], mw[, "Gadus_morhua"])
})

test_that("incompatible water positions are disconnected", {
  tr <- toy_traits()
  mw <- apply_model_metaweb(tr, toy_pars)
  out <- only_correction(mw, tr, water_pos = TRUE)

  # bathydemersal x pelagic-oceanic: removed both ways
  expect_equal(out["Mallotus_villosus", "Lycodes_esmarkii"], 0)
  expect_equal(out["Lycodes_esmarkii", "Mallotus_villosus"], 0)
  # benthopelagic x benthopelagic: kept
  expect_equal(out["Clupea_harengus", "Gadus_morhua"], mw["Clupea_harengus", "Gadus_morhua"])
  expect_true(all(diag(out) == 0))
})

test_that("vectorised water-position correction matches the original pairwise rules", {
  rules <- list(c("bathydemersal", "pelagic-oceanic"), c("bathydemersal", "pelagic-neritic"),
                c("bathydemersal", "reef-associated"), c("bathydemersal", "pelagic"),
                c("bathypelagic", "pelagic-oceanic"), c("bathypelagic", "pelagic-neritic"),
                c("bathypelagic", "reef-associated"), c("pelagic-oceanic", "reef-associated"),
                c("pelagic-oceanic", "demersal"), c("pelagic-oceanic", "benthopelagic"))
  incompatible <- function(a, b) {
    any(vapply(rules, function(p) (a == p[1] && b == p[2]) || (a == p[2] && b == p[1]), logical(1)))
  }
  habitats <- unique(unlist(rules))

  set.seed(1)
  h <- sample(habitats, 30, replace = TRUE)
  m <- matrix(runif(900), 30)
  expected <- m
  for (i in 1:30) for (j in 1:30) if (incompatible(h[i], h[j])) expected[i, j] <- 0

  expect_equal(Distri_correction(h, m), expected)
})

test_that("small fish and producer nodes are handled", {
  out <- toy_metaweb()

  expect_equal(dim(out), c(10, 10))
  expect_identical(tail(rownames(out), 2), c("PrimaryProducer", "SecondaryProducer"))
  # stickleback (5 cm) eats only secondary producers
  expect_equal(sum(out[, "Gasterosteus_aculeatus"]), 1)
  expect_equal(out["SecondaryProducer", "Gasterosteus_aculeatus"], 1)
  expect_equal(unname(out[producer_names, "SecondaryProducer"]), c(1, 1))
  expect_equal(sum(out[, "PrimaryProducer"]), 0)
})

test_that("big low-trophic-level fish eat only secondary producers", {
  tr <- toy_traits()
  tr["Somniosus_microcephalus", "TrophicLevel"] <- 3.2
  out <- toy_metaweb(tr)
  expect_equal(unname(out[, "Somniosus_microcephalus"]), c(rep(0, 9), 1))

  mw <- apply_model_metaweb(tr, toy_pars)
  expect_warning(suppressMessages(correct_metaweb_fish(mw, tr, PP_SP_add = FALSE)), "producer nodes")
})

test_that("traits can contain extra taxa but must cover the metaweb", {
  tr <- toy_traits()
  mw_small <- apply_model_metaweb(tr[1:5, ], toy_pars)
  expect_no_error(suppressMessages(correct_metaweb_fish(mw_small, tr)))

  mw <- apply_model_metaweb(tr, toy_pars)
  expect_error(correct_metaweb_fish(mw, tr[1:5, ]), "not row names")
  expect_error(correct_metaweb_fish(toy_metaweb(), tr), "already contains producer")
})

test_that("missing traits give a warning, not an error", {
  tr <- toy_traits()
  tr["Clupea_harengus", "DemersPelag"] <- NA
  mw <- apply_model_metaweb(tr, toy_pars)
  expect_warning(suppressMessages(correct_metaweb_fish(mw, tr)), "DemersPelag")
})

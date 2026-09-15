test_that("set_herbivores creates the Herbivore column and applies changes", {
  tr <- toy_traits()
  tr["Boreogadus_saida", "TrophicLevel"] <- 2

  out <- set_herbivores(tr)
  expect_identical(rownames(out)[out$Herbivore], "Boreogadus_saida")

  out <- set_herbivores(tr, add = "Clupea_harengus", remove = "Boreogadus_saida")
  expect_identical(rownames(out)[out$Herbivore], "Clupea_harengus")
  # further changes keep earlier ones
  out <- set_herbivores(out, add = "Mallotus_villosus")
  expect_setequal(rownames(out)[out$Herbivore], c("Clupea_harengus", "Mallotus_villosus"))

  expect_error(set_herbivores(tr, add = "Not_a_fish"), "Not_a_fish")
  expect_error(set_herbivores(tr, add = "Gadus_morhua", remove = "Gadus_morhua"), "both")
})

test_that("correct_metaweb_fish lists herbivores and follows the Herbivore column", {
  tr <- toy_traits()
  tr["Boreogadus_saida", "TrophicLevel"] <- 2
  mw <- apply_model_metaweb(tr, toy_pars)

  messages <- paste(capture_messages(correct_metaweb_fish(mw, tr)), collapse = "")
  expect_match(messages, "herbivores .*Boreogadus_saida")
  expect_match(messages, "set_herbivores")

  # removed from herbivores: keeps its fish prey, no primary producer link
  fixed <- set_herbivores(tr, remove = "Boreogadus_saida")
  out <- suppressMessages(correct_metaweb_fish(mw, fixed, water_pos = FALSE, small_fish = FALSE))
  expect_gt(sum(out[toy_taxa, "Boreogadus_saida"]), 0)
  expect_equal(out["PrimaryProducer", "Boreogadus_saida"], 0)

  # added as herbivore: loses its fish prey, eats primary producers
  added <- set_herbivores(tr, add = "Gadus_morhua")
  out <- suppressMessages(correct_metaweb_fish(mw, added))
  expect_equal(sum(out[toy_taxa, "Gadus_morhua"]), 0)
  expect_equal(out["PrimaryProducer", "Gadus_morhua"], 1)
})

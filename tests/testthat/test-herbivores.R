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

test_that("get_herbivores_fishbase uses diet items, then feeding type", {
  local_fake_fishbase()
  res <- get_herbivores_fishbase(c("Chelon_auratus", "Sarpa_salpa", "Gadus_morhua", "Diplodus_sargus",
                                   "Diplodus", "Clupea_harengus"))
  # Chelon: 3/4 plant or detritus records; Sarpa 2/3; Gadus: 3 animal records
  expect_equal(res$plant_share[1:3], c(0.75, 2 / 3, 0))
  expect_identical(res$Herbivore, c(TRUE, TRUE, FALSE, TRUE, TRUE, NA))
  expect_identical(res$Herbivore_source,
                   c("diet_items", "diet_items", "diet_items", "feeding_type", "feeding_type", NA))
  # "variable" is not informative
  expect_true(is.na(get_herbivores_fishbase("Diplodus_vulgaris")$Herbivore))
  # genus pools its species: Diplodus feeding types are a tie, the first alphabetically wins
  expect_identical(res["Diplodus", "feeding_type"], "grazing on aquatic plants")
  expect_identical(get_herbivores_fishbase("Sarpa_salpa", min_plant_share = 0.7)$Herbivore, FALSE)
  expect_identical(get_herbivores_fishbase("Sarpa_salpa", min_records = 4)$Herbivore_source, NA_character_)
})

test_that("set_herbivores(source = 'fishbase') falls back on the trophic level and records sources", {
  local_fake_fishbase()
  tr <- data.frame(TrophicLevel = c(2.4, 2.0, 4.4, 3.1, 2.1),
                   row.names = c("Chelon_auratus", "Sarpa_salpa", "Gadus_morhua", "Diplodus_sargus", "Other_fish"))
  expect_message(out <- set_herbivores(tr, source = "fishbase"), "diet_items: 2")
  expect_identical(out$Herbivore, c(TRUE, TRUE, FALSE, TRUE, TRUE))
  expect_identical(out$Herbivore_source, c("diet_items", "diet_items", "diet_items", "feeding_type", "trophic_level"))

  out <- suppressMessages(set_herbivores(tr, source = "fishbase", remove = "Diplodus_sargus"))
  expect_identical(out["Diplodus_sargus", c("Herbivore", "Herbivore_source")],
                   data.frame(Herbivore = FALSE, Herbivore_source = "manual", row.names = "Diplodus_sargus"))
  # trophic-level rule, as before
  out <- set_herbivores(tr)
  expect_identical(out$Herbivore, c(FALSE, TRUE, FALSE, FALSE, TRUE))
  expect_identical(unique(out$Herbivore_source), "trophic_level")
})

test_that("the correction log records the source of each herbivory change", {
  tr <- toy_traits()
  tr <- set_herbivores(tr, add = "Boreogadus_saida")
  mw <- apply_model_metaweb(tr, toy_pars)
  out <- suppressMessages(correct_metaweb_fish(mw, tr))
  log <- attr(out, "corrections")
  herb <- log[log$correction == "herbivory", ]
  expect_gt(nrow(herb), 0)
  expect_true(all(herb$source == "manual"))
  expect_identical(log$source[log$correction == "water_position"][1], "DemersPelag")
  expect_false(anyNA(log$source))
  expect_true("manual" %in% attr(out, "data_sources"))

  links <- compare_metawebs(mw, out, tr, level = "link", min_proba = 0)
  expect_true(all(links$source[links$correction %in% "herbivory"] == "manual"))
  expect_true("predator_Herbivore_source" %in% names(links))
})

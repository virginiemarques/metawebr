test_that("names are cleaned before matching", {
  x <- c("Diplodus_sp.", "Diplodus spp", "Gobius sp. ABC123", "Gadus cf. morhua", "Salmo trutta trutta",
         "Pomatoschistus minutus/lozanoi", "Pomatoschistus minutus/Gobius niger", "uncultured fish",
         "Salmo salar x Salmo trutta", "Gadus morhua Linnaeus, 1758", "gadus morhua")
  res <- clean_taxon_names(x)
  expect_identical(res$query, c("Diplodus", "Diplodus", "Gobius", "Gadus morhua", "Salmo trutta",
                                "Pomatoschistus", "Pomatoschistus minutus/Gobius niger", "Uncultured fish",
                                "Salmo salar x Salmo trutta", "Gadus morhua", "Gadus morhua"))
  expect_identical(res$not_a_taxon, c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("harmonise_taxa resolves accepted names, synonyms, genera and families with FishBase", {
  local_fake_fishbase()
  taxa <- c("Gadus_morhua", "Liza aurata", "Diplodus sp.", "Sparidae", "Sargus vulgaris",
            "Diplodus sargus", "Nonexistus fish", "uncultured fish", "Liza aurata")
  msg <- capture_messages(res <- harmonise_taxa(taxa, worms = FALSE))
  expect_match(paste(msg, collapse = ""), "Check by hand")

  expect_identical(res$input, unique(taxa))
  expect_identical(res$accepted, c("Gadus_morhua", "Chelon_auratus", "Diplodus", "Sparidae", NA,
                                   "Diplodus_sargus", NA, NA))
  expect_identical(res$status, c("accepted", "synonym", "accepted", "accepted", "ambiguous",
                                 "accepted", "unresolved", "not_a_taxon"))
  expect_identical(res$rank, c("species", "species", "genus", "family", NA, "species", NA, NA))
  expect_identical(res$SpecCode, c(1L, 2L, NA, NA, NA, 4L, NA, NA))
  # misapplied names are ignored: Diplodus sargus stays itself
  expect_match(res$note[5], "Diplodus sargus, Diplodus vulgaris")
})

test_that("user lookup overrides FishBase", {
  local_fake_fishbase()
  lookup <- data.frame(input = c("Liza_aurata", "Nonexistus fish"), accepted = c("Chelon auratus", NA))
  res <- suppressMessages(harmonise_taxa(c("Liza aurata", "Nonexistus fish"), lookup = lookup, worms = FALSE))
  expect_identical(res$accepted, c("Chelon_auratus", NA))
  expect_identical(res$status, c("user", "not_a_taxon"))
  expect_error(harmonise_taxa("x", lookup = data.frame(a = 1), worms = FALSE), "lookup")
})

test_that("WoRMS resolves names unknown to FishBase, and spelling mistakes", {
  skip_if_not_installed("worrms")
  local_fake_fishbase()
  local_mocked_bindings(worms_match = function(names) {
    known <- data.frame(q = c("Gadus morua", "Phoca vitulina", "Sargus vulgaris"),
                        worms_name = c("Gadus morhua", "Phoca vitulina", "Diplodus vulgaris"),
                        AphiaID = c(126436L, 137084L, 127054L), match_type = c("phonetic", "exact", "exact"))
    known[match(names, known$q), c("worms_name", "AphiaID", "match_type")]
  })
  res <- suppressMessages(harmonise_taxa(c("Gadus morua", "Phoca vitulina", "Sargus vulgaris", "Abc def")))
  expect_identical(res$accepted, c("Gadus_morhua", NA, "Diplodus_vulgaris", NA))
  expect_identical(res$status, c("worms_fuzzy", "not_in_fishbase", "worms_exact", "unresolved"))
  expect_identical(res$source, c("worms", NA, "worms", NA))
  expect_identical(res$AphiaID, c(126436L, 137084L, 127054L, NA))
})

test_that("rename_taxa replaces names and keeps the originals", {
  det <- data.frame(taxon = c("Liza aurata", "Chelon_auratus", "Gobius sp.", "uncultured fish"),
                    reads = 1:4)
  harm <- data.frame(input = det$taxon, accepted = c("Chelon_auratus", "Chelon_auratus", "Gobius", NA))
  expect_message(out <- rename_taxa(det, "taxon", harm), "1 accepted names gather")
  expect_identical(out$taxon, c("Chelon_auratus", "Chelon_auratus", "Gobius", "uncultured_fish"))
  expect_identical(out$taxon_original, det$taxon)

  out <- suppressMessages(rename_taxa(det, "taxon", harm, drop_unresolved = TRUE))
  expect_identical(nrow(out), 3L)
  expect_warning(suppressMessages(rename_taxa(det, "taxon", harm[1:3, ])), "uncultured fish")
})

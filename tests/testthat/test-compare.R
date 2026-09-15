herbivore_toy <- function() {
  tr <- toy_traits()
  tr["Boreogadus_saida", "TrophicLevel"] <- 2
  mw <- apply_model_metaweb(tr, toy_pars)
  list(traits = tr, mw = mw, mw_corr = suppressMessages(correct_metaweb_fish(mw, tr)))
}

# Metaweb `m` extended with zero rows/columns to the taxa of `template`
expand_like <- function(m, template) {
  out <- matrix(0, nrow(template), ncol(template), dimnames = dimnames(template))
  out[rownames(m), colnames(m)] <- m
  out
}

test_that("pair table lists each changed species pair once, with both directions", {
  d <- herbivore_toy()
  pairs <- compare_metawebs(d$mw, d$mw_corr, min_proba = 0)

  unordered <- ifelse(pairs$taxon_1 < pairs$taxon_2,
                      paste(pairs$taxon_1, pairs$taxon_2), paste(pairs$taxon_2, pairs$taxon_1))
  expect_false(anyDuplicated(unordered) > 0)

  changed <- abs(unclass(d$mw_corr) - expand_like(d$mw, d$mw_corr)) > 1e-12
  expect_equal(nrow(pairs), sum((changed | t(changed))[upper.tri(changed, diag = TRUE)]))

  # both directions of an incompatible habitat pair are in one row
  row <- pairs[(pairs$taxon_1 == "Mallotus_villosus" & pairs$taxon_2 == "Lycodes_esmarkii") |
               (pairs$taxon_1 == "Lycodes_esmarkii" & pairs$taxon_2 == "Mallotus_villosus"), ]
  expect_equal(nrow(row), 1)
  expect_equal(c(row$t1_eats_t2_correction, row$t2_eats_t1_correction), rep("water_position", 2))
  expect_equal(row$t1_eats_t2_before, d$mw[row$taxon_2, row$taxon_1])
  expect_equal(row$t2_eats_t1_before, d$mw[row$taxon_1, row$taxon_2])
  expect_equal(c(row$t1_eats_t2_after, row$t2_eats_t1_after), c(0, 0))

  # cannibalism: a single direction
  self <- pairs[pairs$taxon_1 == "Gadus_morhua" & pairs$taxon_2 == "Gadus_morhua", ]
  expect_equal(self$t1_eats_t2_correction, "cannibalism")
  expect_true(is.na(self$t2_eats_t1_before))
})

test_that("pair table is filtered, sorted and has traits", {
  d <- herbivore_toy()
  pairs <- compare_metawebs(d$mw, d$mw_corr, d$traits)

  expect_true(all(c("t1_TrophicLevel", "t2_DemersPelag") %in% names(pairs)))
  relevant <- function(b, a, corr) !is.na(b) & abs(a - b) > 1e-12 & pmax(b, a) >= 0.01
  expect_true(all(relevant(pairs$t1_eats_t2_before, pairs$t1_eats_t2_after) |
                  relevant(pairs$t2_eats_t1_before, pairs$t2_eats_t1_after)))
  expect_lt(nrow(pairs), nrow(compare_metawebs(d$mw, d$mw_corr, min_proba = 0)))
})

test_that("link table lists every corrected link with the correction responsible", {
  d <- herbivore_toy()
  links <- compare_metawebs(d$mw, d$mw_corr, level = "link", min_proba = 0)

  changed <- abs(unclass(d$mw_corr) - expand_like(d$mw, d$mw_corr)) > 1e-12
  expect_equal(nrow(links), sum(changed))
  expect_false(anyNA(links$correction))
  expect_setequal(links$correction,
                  c("herbivory", "water_position", "cannibalism", "small_fish", "producers"))

  herb <- links[links$predator == "Boreogadus_saida" & links$prey %in% toy_taxa, ]
  expect_true(all(herb$correction == "herbivory" & herb$change == "removed"))
  expect_equal(herb$before, unname(d$mw[herb$prey, "Boreogadus_saida"]))

  added <- links[links$change == "added", ]
  expect_true(all(added$correction == "producers" & added$prey %in% producer_names))

  filtered <- compare_metawebs(d$mw, d$mw_corr, d$traits, level = "link")
  expect_true(all(pmax(filtered$before, filtered$after) >= 0.01))
  expect_true(all(diff(filtered$before) <= 0))
  expect_equal(unique(filtered$predator_TrophicLevel[filtered$predator == "Boreogadus_saida"]), 2)
})

test_that("changes not made by correct_metaweb_fish get no correction", {
  d <- herbivore_toy()

  bin <- correct_metaweb_proba(d$mw_corr, threshold = 0.5)
  expect_null(attr(bin, "corrections"))
  links <- compare_metawebs(d$mw_corr, bin, level = "link", min_proba = 0)
  expect_gt(nrow(links), 0)
  expect_true(all(is.na(links$correction)))

  # a link edited by hand before the comparison is no longer attributed
  edited <- d$mw
  edited["Mallotus_villosus", "Lycodes_esmarkii"] <- 0.123
  links <- compare_metawebs(edited, d$mw_corr, level = "link", min_proba = 0)
  hand <- links$prey == "Mallotus_villosus" & links$predator == "Lycodes_esmarkii"
  expect_true(is.na(links$correction[hand]))
  expect_false(all(is.na(links$correction)))
})

test_that("taxon-level summary adds up", {
  d <- herbivore_toy()
  s <- compare_metawebs(d$mw, d$mw_corr, d$traits, level = "taxon")

  expect_identical(s$taxon, rownames(d$mw_corr))
  expect_equal(s$prey_before[1:8], unname(colSums(d$mw)))
  expect_equal(s$prey_after, unname(colSums(d$mw_corr)))
  change_cols <- grep("^prey_change_", names(s), value = TRUE)
  expect_equal(s$prey_after - s$prey_before, rowSums(s[change_cols]))

  bor <- s[s$taxon == "Boreogadus_saida", ]
  expect_equal(bor$prey_change_herbivory, -sum(d$mw[, "Boreogadus_saida"]))
  expect_equal(bor$TrophicLevel, 2)
  expect_error(compare_metawebs(d$mw, d$mw_corr, level = "species"), "should be one of")
})

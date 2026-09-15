fake_fishbase <- function() {
  data.frame(
    Species = c("Gadus_morhua", "Gadus_chalcogrammus", "Clupea_harengus", "Clupea_pallasii"),
    Genus = c("Gadus", "Gadus", "Clupea", "Clupea"),
    Family = c("Gadidae", "Gadidae", "Clupeidae", "Clupeidae"),
    CommonLengthEstim = c(100, 50, 30, NA),
    TrophicLevel = c(4.1, 3.5, 3.4, 3.2),
    DemersPelag = c("benthopelagic", "benthopelagic", "pelagic-neritic", NA)
  )
}

test_that("get_traits works with any species column name (regression)", {
  local_mocked_bindings(fishbase_traits_table = function() fake_fishbase())
  presence <- data.frame(Species = c("Gadus morhua", "Clupea_harengus", "Gadus morhua", "Nonexistus_fish"))

  expect_warning(res <- get_traits(presence, column_species = "Species"), "Nonexistus_fish")
  expect_identical(rownames(res), c("Gadus_morhua", "Clupea_harengus", "Nonexistus_fish"))
  expect_named(res, c("Species", "CommonLengthEstim", "TrophicLevel", "DemersPelag"))
  expect_equal(res["Gadus_morhua", "CommonLengthEstim"], 100)
})

test_that("infer_traits averages traits at genus and family level", {
  taxa <- c("Gadus_morhua", "Gadus", "Clupeidae", "Nonexistus")
  expect_warning(res <- infer_traits(taxa, fishbase_table = fake_fishbase()), "Nonexistus")

  expect_identical(rownames(res), taxa)
  expect_equal(res$rank, c("species", "genus", "family", NA))
  expect_equal(res$n_species, c(1, 2, 2, 0))
  expect_equal(res$CommonLengthEstim, c(100, 75, 30, NA))
  expect_equal(res$TrophicLevel, c(4.1, 3.8, 3.3, NA))
  expect_equal(res$DemersPelag, c("benthopelagic", "benthopelagic", "pelagic-neritic", NA))
})

test_that("clean_traits fills missing values by genus and drops incomplete taxa", {
  tr <- data.frame(
    CommonLengthEstim = c(100, NA, 30, NA),
    TrophicLevel = c(4, 3, NA, 3),
    DemersPelag = c("demersal", NA, "pelagic", "pelagic"),
    row.names = c("Gadus_morhua", "Gadus_ogac", "Clupea_harengus", "Mallotus_villosus")
  )
  out <- suppressMessages(clean_traits(tr))

  expect_equal(out["Gadus_ogac", "CommonLengthEstim"], 100)
  expect_equal(out["Gadus_ogac", "DemersPelag"], "demersal")
  expect_identical(rownames(out), c("Gadus_morhua", "Gadus_ogac"))

  kept <- suppressMessages(clean_traits(tr, fill_by_genus = FALSE, remove_incomplete = FALSE))
  expect_identical(kept, tr)
})

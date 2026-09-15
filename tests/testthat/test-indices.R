test_that("top, basal and intermediate fractions are computed (regression: always NA)", {
  nm <- c("A", "B", "C", "D", "PrimaryProducer", "SecondaryProducer")
  w <- matrix(0, 6, 6, dimnames = list(nm, nm))
  w["PrimaryProducer", "A"] <- 1
  w["SecondaryProducer", "B"] <- 1
  w["PrimaryProducer", "SecondaryProducer"] <- 1
  w["A", "C"] <- 1
  w["B", "C"] <- 1
  w["C", "D"] <- 1

  ind <- Calc_indic_proba(w)
  # without producers: prey counts A0 B0 C2 D1, predator counts A1 B1 C1 D0
  expect_equal(unname(ind[c("Ntop", "Nbas", "Nint", "Vul", "Gen")]), c(0.25, 0.75, 0, 1, 1.5))
  expect_equal(unname(ind[c("Species", "Link")]), c(6, 6))
  # trophic levels: PP 1, SP 2, A 2, B 3, C 3.5, D 4.5
  expect_equal(unname(ind[c("TL_moy", "length_chain")]), c(16 / 6, 5))
})

test_that("path statistics always have five length classes (regression)", {
  nm <- paste0("s", 1:9)
  chain <- matrix(0, 9, 9, dimnames = list(nm, nm))
  for (i in 1:8) chain[i, i + 1] <- 1  # paths up to 8 links long

  long <- Calc_indic_proba(chain)
  short <- Calc_indic_proba(chain[1:3, 1:3])
  expect_identical(names(long), names(short))
  expect_false(anyNA(long[paste0("Shortest_path_", 1:5)]))
  expect_true(all(is.na(short[paste0("Shortest_path_", 3:5)])))
})

test_that("get_indic_cells returns one row per site, including single-taxon sites", {
  bin <- correct_metaweb_proba(toy_metaweb(), threshold = 0.5)
  det <- data.frame(station = c("s1", "s1", "s1", "s2", "s2", "s3"),
                    taxon = c("Gadus_morhua", "Clupea_harengus", "Mallotus_villosus",
                              "Boreogadus_saida", "Gadus_morhua", "Clupea_harengus"),
                    count = 1)
  pa <- make_species_station_matrix(det, "station", "taxon", "count")

  res <- get_indic_cells(pa, bin)
  expect_s3_class(res, "data.frame")
  expect_identical(res$site, c("s1", "s2", "s3"))
  expect_equal(res$Species, c(5, 4, 3))
})

test_that("taxa absent from the metaweb are dropped with a warning", {
  bin <- correct_metaweb_proba(toy_metaweb(), threshold = 0.5)
  pa <- toy_sites()
  pa_extra <- cbind(pa, Unknown_fish = c(1, 0, 0))

  expect_warning(res <- get_indic_cells(pa_extra, bin), "Unknown_fish")
  expect_equal(res, get_indic_cells(pa, bin))
})

test_that("make_species_station_matrix builds a binary matrix with producers", {
  det <- data.frame(station = c("s1", "s1", "s2"), taxon = c("A", "B", "A"), reads = c(10, 0, 3))
  m <- make_species_station_matrix(det, "station", "taxon", "reads")
  expected <- matrix(c(1, 1, 0, 0, 1, 1, 1, 1), 2, 4,
                     dimnames = list(c("s1", "s2"), c("A", "B", "PrimaryProducer", "SecondaryProducer")))
  expect_equal(m, expected)
  expect_error(make_species_station_matrix(det, "station", "taxon", "count"), "must contain columns")
})

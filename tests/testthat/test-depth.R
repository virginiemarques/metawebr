test_that("add_depth_range reads FishBase and df_depth_env, keeping existing values", {
  local_fake_fishbase()
  tr <- data.frame(TrophicLevel = c(4.4, 3.1, 3), row.names = c("Gadus_morhua", "Diplodus", "Nonexistus_fish"))
  expect_message(out <- add_depth_range(tr), "for 2 taxa; 1 taxa without")
  expect_equal(out$DepthMin, c(150, 0, NA))
  expect_equal(out$DepthMax, c(600, 160, NA))  # genus: widest range of its species
  expect_identical(out$Depth_source, c("fishbase", "fishbase", NA))

  sp <- gsub(" ", "_", df_depth_env$Species[1])
  tr2 <- data.frame(TrophicLevel = c(4.4, 3), row.names = c("Gadus_morhua", sp))
  out2 <- suppressMessages(add_depth_range(add_depth_range(tr2), "df_depth_env"))
  expect_identical(out2$Depth_source, c("fishbase", "df_depth_env"))
  expect_equal(out2[sp, "DepthMax"], df_depth_env$Depth_max[1])
})

test_that("depth_overlap removes links between taxa that never meet, and logs them", {
  tr <- toy_traits()
  tr$DepthMin <- c(0, 0, 0, 500, 0, 0, 0, NA)
  tr$DepthMax <- c(200, 200, 200, 1000, 200, 2000, 50, NA)
  tr$Depth_source <- "test"
  mw <- apply_model_metaweb(tr, toy_pars)

  base <- suppressMessages(correct_metaweb_fish(mw, tr))
  out <- suppressMessages(suppressWarnings(correct_metaweb_fish(mw, tr, depth_overlap = TRUE)))
  deep <- "Reinhardtius_hippoglossoides"
  shallow <- c("Gadus_morhua", "Clupea_harengus", "Mallotus_villosus", "Boreogadus_saida",
               "Gasterosteus_aculeatus")
  expect_true(all(out[shallow, deep] == 0))
  expect_true(all(out[deep, shallow] == 0))
  # overlapping or unknown ranges are untouched
  expect_equal(out["Clupea_harengus", "Gadus_morhua"], base["Clupea_harengus", "Gadus_morhua"])
  expect_equal(out[, "Lycodes_esmarkii"], base[, "Lycodes_esmarkii"])

  expect_warning(suppressMessages(correct_metaweb_fish(mw, tr, depth_overlap = TRUE)), "missing `DepthMin_DepthMax`.*Lycodes")
  log <- attr(out, "corrections")
  expect_true(all(log$source[log$correction == "depth"] == "DepthMin/DepthMax (test)"))

  # a tolerance bridges the gaps (300 m, and 450 m to Gasterosteus)
  tol <- suppressMessages(suppressWarnings(correct_metaweb_fish(mw, tr, depth_overlap = TRUE, depth_tolerance = 300)))
  expect_equal(tol["Gasterosteus_aculeatus", deep], 0)
  others <- setdiff(shallow, "Gasterosteus_aculeatus")
  expect_equal(tol[others, deep], base[others, deep])
  tol <- suppressMessages(suppressWarnings(correct_metaweb_fish(mw, tr, depth_overlap = TRUE, depth_tolerance = 450)))
  expect_equal(tol, base, ignore_attr = TRUE)
  expect_error(correct_metaweb_fish(mw, toy_traits(), depth_overlap = TRUE), "DepthMin")
})

test_that("depth_apart compares ranges in both directions", {
  apart <- depth_apart(c(0, 100, 0), c(50, 200, NA))
  expect_identical(apart, matrix(c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 3))
})

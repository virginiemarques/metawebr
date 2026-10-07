test_that("fb_table downloads once per session and can use a disk cache", {
  clear_fishbase_cache()
  withr_dir <- tempfile("fbcache")
  old <- options(metawebr.cache_dir = withr_dir)
  on.exit({ options(old); clear_fishbase_cache(); unlink(withr_dir, recursive = TRUE) })

  calls <- 0
  local_mocked_bindings(fishbase_download = function(table) {
    calls <<- calls + 1
    data.frame(SpecCode = 1:3)
  })
  x <- fb_table("species")
  y <- fb_table("species")
  expect_equal(calls, 1)
  expect_identical(x, y)
  expect_true(file.exists(file.path(withr_dir, "fishbase_species.rds")))

  info <- fishbase_cache_info()
  expect_identical(info$table, "species")
  expect_identical(info$location, "memory")

  # new session: read from disk, no download
  clear_fishbase_cache()
  expect_identical(fishbase_cache_info()$location, "disk")
  z <- fb_table("species")
  expect_equal(calls, 1)
  expect_equal(z, x)

  clear_fishbase_cache(disk = TRUE)
  expect_false(file.exists(file.path(withr_dir, "fishbase_species.rds")))
  expect_error(fb_table("nope"))
})

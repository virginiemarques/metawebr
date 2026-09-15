toy_tree_traits <- function() {
  set.seed(10)
  n <- 30
  tree <- ape::rcoal(n, tip.label = paste0("Genus", rep(1:6, each = 5), "_sp", 1:n))
  # traits with a phylogenetic signal: close relatives have similar values
  tl <- 2.5 + as.numeric(ape::rTraitCont(tree, sigma = 0.5))
  traits <- data.frame(
    CommonLengthEstim = exp(3 + as.numeric(ape::rTraitCont(tree, sigma = 0.5))),
    TrophicLevel = tl,
    DemersPelag = ifelse(tl > median(tl), "demersal", "reef-associated"),
    row.names = tree$tip.label
  )
  list(tree = tree, traits = traits)
}

test_that("phylogenetic eigenvectors match PVR::PVRdecomp()", {
  skip_if_not_installed("ape")
  skip_if_not_installed("PVR")
  tree <- toy_tree_traits()$tree
  expect_equal(unname(phylo_eigenvectors(tree, 5)),
               unname(PVR::PVRdecomp(tree)@Eigen$vectors[, 1:5]), tolerance = 1e-8)
})

test_that("impute_traits_phylo fills only missing values of species in the tree", {
  skip_if_not_installed("ape")
  skip_if_not_installed("missForest")
  d <- toy_tree_traits()
  traits <- d$traits
  traits[c("Genus1_sp1", "Genus2_sp7"), "CommonLengthEstim"] <- NA
  traits["Genus3_sp12", c("TrophicLevel", "DemersPelag")] <- NA
  traits["Outside_tree", ] <- list(NA, 3, "demersal")
  trees <- list(d$tree, d$tree)

  set.seed(1)
  expect_warning(
    res <- suppressMessages(impute_traits_phylo(traits, trees = trees, n_eigen = 5, n_runs = 2,
                                                maxiter = 3, ntree = 20)),
    "Outside_tree"
  )

  observed <- !is.na(traits$CommonLengthEstim)
  expect_equal(res$CommonLengthEstim[observed], traits$CommonLengthEstim[observed])
  expect_false(anyNA(res[rownames(traits) != "Outside_tree", c("CommonLengthEstim", "TrophicLevel", "DemersPelag")]))
  expect_true(is.na(res["Outside_tree", "CommonLengthEstim"]))
  expect_true(res["Genus3_sp12", "DemersPelag"] %in% c("demersal", "reef-associated"))
  expect_equal(res["Genus3_sp12", "imputed"], "TrophicLevel, DemersPelag")
  expect_equal(res["Genus1_sp1", "imputed"], "CommonLengthEstim")
  expect_true(is.na(res["Genus4_sp16", "imputed"]))

  set.seed(1)
  again <- suppressWarnings(suppressMessages(
    impute_traits_phylo(traits, trees = trees, n_eigen = 5, n_runs = 2, maxiter = 3, ntree = 20)
  ))
  expect_identical(res, again)
})

test_that("infer_traits(method = '100matrix') imputes after the FishBase lookup", {
  skip_if_not_installed("ape")
  skip_if_not_installed("missForest")
  d <- toy_tree_traits()
  fishbase <- data.frame(Species = rownames(d$traits), Genus = sub("_.*", "", rownames(d$traits)),
                         Family = "Fakeidae", d$traits)
  fishbase$CommonLengthEstim[1] <- NA

  set.seed(2)
  res <- suppressMessages(infer_traits(rownames(d$traits), fishbase_table = fishbase, method = "100matrix",
                                       trees = d$tree, n_eigen = 5, n_runs = 2, maxiter = 3, ntree = 20))
  expect_false(anyNA(res$CommonLengthEstim))
  expect_equal(res$imputed[1], "CommonLengthEstim")
  expect_error(infer_traits("Genus1_sp1", fishbase_table = fishbase, method = "other"), "should be one of")
})

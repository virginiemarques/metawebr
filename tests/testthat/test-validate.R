web3 <- function() {
  matrix(c(0, 0.9, 0.2,
           0, 0, 0.8,
           0, 0, 0), 3, byrow = TRUE, dimnames = rep(list(c("A", "B", "C")), 2))
}

test_that("auc and confusion scores are correct", {
  expect_equal(auc(c(0.9, 0.8, 0.1, 0.2), c(1, 1, 0, 0)), 1)
  expect_equal(auc(c(0.1, 0.2, 0.9, 0.8), c(1, 1, 0, 0)), 0)
  expect_equal(auc(c(0.5, 0.5), c(1, 0)), 0.5)
  expect_true(is.na(auc(1:3, c(1, 1, 1))))
  s <- confusion_scores(c(TRUE, FALSE, TRUE, FALSE), c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(unname(s), c(0.5, 0.5, 0, 0.5))
})

test_that("validate_metaweb scores observed links against pseudo-absences", {
  obs <- data.frame(predator = c("B", "C"), prey = c("A", "B"))
  res <- suppressMessages(validate_metaweb(web3(), obs))
  # pairs: B eats A (0.9, obs), C (0); C eats A (0.2), B (0.8, obs)
  expect_equal(res$summary$n_pairs, 4)
  expect_equal(res$summary$n_links, 2)
  expect_equal(res$summary$AUC, 1)
  expect_equal(res$summary$TSS, 1)
  expect_identical(res$per_predator$missed, c("", ""))

  res <- suppressMessages(validate_metaweb(web3(), obs, threshold = 0.85))
  expect_equal(res$summary$sensitivity, 0.5)
  expect_identical(res$per_predator$missed[res$per_predator$predator == "C"], "B")

  # explicit absences: only the listed pairs
  obs01 <- data.frame(predator = c("B", "C"), prey = c("A", "A"), interaction = c(1, 0))
  expect_equal(suppressMessages(validate_metaweb(web3(), obs01))$summary$n_pairs, 2)

  # binary metaweb: no AUC
  bin <- (web3() >= 0.5) * 1
  expect_true(is.na(suppressMessages(validate_metaweb(bin, obs))$summary$AUC))
  expect_error(suppressMessages(validate_metaweb(web3(), data.frame(predator = "X", prey = "Y"))),
               "No observed link")
})

test_that("validate_metaweb warns when validation data were used to correct the metaweb", {
  mw <- web3()
  attr(mw, "data_sources") <- c("diet_items", "DemersPelag")
  obs <- data.frame(predator = "B", prey = "A")
  attr(obs, "source") <- "fishbase_fooditems"
  expect_warning(validate_metaweb(mw, obs), "not independent")
})

test_that("get_proba can choose the threshold from observed links", {
  obs <- data.frame(predator = c("B", "C"), prey = c("A", "B"))
  th <- suppressMessages(get_proba(web3(), plot_diag = FALSE, method = "tss", links = obs))
  expect_true(th > 0.2 && th <= 0.8)
  bin <- suppressMessages(correct_metaweb_proba(web3(), threshold = NULL, plot_diag = FALSE,
                                                method = "tss", links = obs))
  expect_equal(sum(bin), 2)
  expect_error(get_proba(web3(), method = "tss", plot_diag = FALSE), "links")
  expect_error(get_proba(web3(), plot_diag = FALSE), "df_traits")
})

test_that("get_diet_links maps FishBase prey to species and genera", {
  local_fake_fishbase()
  expect_message(links <- get_diet_links(c("Gadus_morhua", "Clupea_harengus", "Diplodus", "PrimaryProducer")),
                 "2 links")
  expect_identical(links$predator, c("Gadus_morhua", "Gadus_morhua"))
  expect_identical(links$prey, c("Clupea_harengus", "Diplodus"))
  expect_identical(attr(links, "source"), "fishbase_fooditems")
})

test_that("cv_size_model holds out whole studies", {
  set.seed(1)
  # two studies drawn from the same size relationship
  sim <- function(ref, n) {
    pred_len <- 10^runif(n, 1, 2)
    # larger predators eat larger prey (prey k is 10^(k/4) cm)
    k <- pmin(8, pmax(1, round(4 * log10(pred_len) - 3 + sample(-1:1, n, replace = TRUE))))
    prey <- paste0("prey", k)
    prey_len <- 10^(k / 4)
    data.frame(reference = ref, predator = paste0(ref, "_pred", ceiling(log10(pred_len) * 3)),
               prey = prey, standardised_predator_length = pred_len, si_prey_length = prey_len)
  }
  d <- rbind(sim("s1", 60), sim("s2", 60), sim("s3", 2))
  cv <- suppressWarnings(cv_size_model(d, max_time = 1, min_links = 5, prey_pool = "study"))
  expect_identical(cv$folds$group, c("s1", "s2", "s3"))
  expect_true(all(!is.na(cv$folds$AUC[1:2])))
  expect_true(is.na(cv$folds$AUC[3]))  # too few links to be scored
  expect_equal(cv$folds$n_train, c(62, 62, 120))
  expect_identical(cv$overall$model, c("out_of_sample", "in_sample"))
  expect_setequal(unique(cv$pairs$group), c("s1", "s2"))
  expect_error(cv_size_model(d, group = "nope"), "no column")

  # all prey of the data: single-predator studies can be scored
  one <- rbind(d, data.frame(reference = "s4", predator = "solo", prey = c("prey1", "prey2", "prey3"),
                             standardised_predator_length = 20, si_prey_length = 10^(1:3 / 4)))
  cv_all <- suppressWarnings(cv_size_model(one, max_time = 1, min_links = 3, in_sample = FALSE))
  expect_false(is.na(cv_all$folds$AUC[cv_all$folds$group == "s4"]))
  expect_equal(cv_all$folds$n_prey[cv_all$folds$group == "s4"], length(unique(one$prey)))
  expect_true(all(is.na(cv_all$pairs$proba_in_sample)))
})

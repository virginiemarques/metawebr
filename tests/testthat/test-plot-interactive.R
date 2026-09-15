toy_web_bin <- function() {
  (toy_metaweb() >= 0.5) * 1
}

test_that("plot_network_interactive() draws every taxon and link of the web", {
  skip_if_not_installed("visNetwork")
  mw <- toy_web_bin()
  web <- plot_network_interactive(mw)

  expect_s3_class(web, "visNetwork")
  expect_setequal(web$x$nodes$id, rownames(mw))
  expect_equal(nrow(web$x$edges), sum(mw > 0))
  # links go from prey (rows) to predators (columns)
  expect_true(all(mw[cbind(web$x$edges$from, web$x$edges$to)] == 1))
})

test_that("plot_network_interactive() keeps sp_to_keep and the producers", {
  skip_if_not_installed("visNetwork")
  mw <- toy_web_bin()
  keep <- toy_taxa[1:4]
  web <- plot_network_interactive(mw, sp_to_keep = keep)

  expect_setequal(web$x$nodes$id, intersect(c(keep, producer_names), rownames(mw)))
  expect_true(all(c(web$x$edges$from, web$x$edges$to) %in% web$x$nodes$id))
})

test_that("plot_network_interactive() uses the layout of plot_tree_network()", {
  skip_if_not_installed("visNetwork")
  mw <- toy_web_bin()
  set.seed(42)
  web <- plot_network_interactive(mw)
  set.seed(42)
  p <- plot_tree_network(mw)

  expect_equal(order(web$x$nodes$x), order(p$data$x))
  expect_equal(order(-web$x$nodes$y), order(p$data$y))
})

test_that("plot_network_interactive() checks its arguments", {
  skip_if_not_installed("visNetwork")
  mw <- toy_web_bin()
  expect_error(plot_network_interactive(mw, links = "neighbours"), "should be one of")
  expect_error(plot_network_interactive(mw, sp_to_keep = "Not_a_fish"), "None of `sp_to_keep`")
  expect_error(plot_network_interactive(matrix(1, 2, 3)), "must be square")
})

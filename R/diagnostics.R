#' Plot the fit of the niche model for one predator
#'
#' Left panel: size distribution of the observed prey of `target_sp`, of all
#' prey, and the predicted niche (scaled) for the mean size of `target_sp`.
#' Right panel: predicted interaction probability over the whole predator-prey
#' size space, with all observations (those of `target_sp` highlighted). Sizes
#' are on a log10 scale, with axes labelled in cm.
#'
#' @param data Interaction data with the columns `predator`,
#'   `standardised_predator_length` and `si_prey_length`, e.g.
#'   [df_interaction_fish].
#' @param target_sp Predator name, as in `data$predator`.
#' @param pars Model parameters, in any form accepted by
#'   [apply_model_metaweb()].
#'
#' @return `NULL`, invisibly. Called for its plot.
#'
#' @export
get_fig_eval <- function(data, target_sp = "Zeus faber", pars) {
  cols <- c("predator", "standardised_predator_length", "si_prey_length")
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("`data` is missing column(s): %s.", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (!target_sp %in% data$predator) {
    stop(sprintf("`target_sp` '%s' is not a predator in `data`.", target_sp), call. = FALSE)
  }
  pars <- as_pars(pars)

  col_target <- "#eb6834"  # observations of target_sp
  col_all <- "#c3c2b7"     # all prey
  col_model <- "#1c5cab"   # model prediction

  MPrey <- log10(data$si_prey_length)
  MPred <- log10(data$standardised_predator_length)
  is_target <- data$predator == target_sp

  # Share of prey per size bin
  bars_all <- hist_prop(MPrey, 20, "All prey")
  bars_target <- hist_prop(MPrey[is_target], 10, "Observed prey")
  # Outline of the all-prey histogram, drawn over the target bars
  outline_all <- data.frame(size = rep(c(bars_all$xmin, utils::tail(bars_all$xmax, 1)), each = 2),
                            prop = c(0, rep(bars_all$prop, each = 2), 0))

  # Predicted niche for the mean size of the target predator
  target_size <- mean(MPred[is_target])
  seqM <- seq(min(MPrey) - 0.5, max(MPrey), 0.01)
  niche <- data.frame(size = seqM, prop = 0.2 * pLMFitted(seqM, target_size, pars))

  p_niche <- ggplot() +
    geom_rect(data = rbind(bars_all, bars_target),
              aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = .data$prop, fill = .data$group),
              colour = "white", linewidth = 0.3) +
    geom_path(data = outline_all, aes(.data$size, .data$prop), colour = "#898781", linewidth = 0.35) +
    geom_hline(yintercept = 0, colour = "#898781", linewidth = 0.4) +
    geom_line(data = niche, aes(.data$size, .data$prop, colour = "Predicted niche (scaled)"),
              linewidth = 1) +
    scale_fill_manual(values = c(`Observed prey` = col_target, `All prey` = col_all),
                      breaks = c("Observed prey", "All prey"), name = NULL) +
    scale_colour_manual(values = col_model, name = NULL) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
    scale_x_continuous(breaks = log10(size_breaks_cm), labels = size_breaks_cm, minor_breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(xlim = range(MPrey), ylim = c(0, max(0.4, bars_all$prop, bars_target$prop))) +
    labs(title = "Prey size distribution", x = "Prey body size (cm)", y = "Proportion of prey") +
    theme_fig_eval() +
    theme(panel.grid.major.x = element_blank())

  # Interaction probability over the size space
  seqX <- seq(min(MPred), max(MPred), 0.01)
  seqY <- seq(min(MPrey), max(MPrey), 0.01)
  surface <- expand.grid(pred = seqX, prey = seqY)
  surface$proba <- pLMFitted(surface$prey, surface$pred, pars)

  obs <- data.frame(pred = MPred, prey = MPrey,
                    group = ifelse(is_target, "Observed prey", "All observations"))
  obs_all <- obs[!is_target, ]
  obs_target <- obs[is_target, ]

  p_space <- ggplot(mapping = aes(.data$pred, .data$prey)) +
    geom_raster(data = surface, aes(fill = .data$proba)) +
    geom_point(data = obs_all, aes(colour = .data$group), size = 0.6, stroke = 0, alpha = 0.7) +
    # White halo keeps the target points distinct from the other observations
    geom_point(data = obs_target, colour = "white", size = 2.6, stroke = 0) +
    geom_point(data = obs_target, aes(colour = .data$group), size = 1.8, stroke = 0) +
    scale_fill_gradientn(colours = c("#f6f9fd", "#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab"),
                         limits = c(0, 1), name = "Link\nprobability") +
    scale_colour_manual(values = c(`Observed prey` = col_target, `All observations` = "#1f1f1d"),
                        breaks = c("Observed prey", "All observations"), name = NULL) +
    scale_x_continuous(breaks = log10(size_breaks_cm), labels = size_breaks_cm) +
    scale_y_continuous(breaks = log10(size_breaks_cm), labels = size_breaks_cm) +
    coord_cartesian(expand = FALSE) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2.5, alpha = 1)),
           fill = guide_colourbar(order = 2)) +
    labs(title = "Predicted link probability", x = "Predator body size (cm)", y = "Prey body size (cm)") +
    theme_fig_eval() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.justification = "top",
          legend.key.height = grid::unit(1.2, "lines"))

  # Both panels side by side, with aligned plotting areas, under a shared title
  panels <- cbind(ggplotGrob(p_niche), ggplotGrob(p_space), size = "max")
  title <- grid::textGrob(bquote(bold("Predator:") ~ bolditalic(.(target_sp))),
                          x = grid::unit(8, "pt"), hjust = 0,
                          gp = grid::gpar(fontsize = 14, col = "#0b0b0b"))

  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  layout <- grid::grid.layout(2, 1, heights = grid::unit(c(2.4, 1), c("lines", "null")))
  grid::pushViewport(grid::viewport(layout = layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1))
  grid::grid.draw(title)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2))
  grid::grid.draw(panels)
  grid::upViewport(2)

  invisible(NULL)
}

# Axis breaks (cm) for body sizes plotted on a log10 scale
size_breaks_cm <- c(0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000)

# Histogram of `x` as the share of values per bin, as a data frame of bars
hist_prop <- function(x, breaks, group) {
  h <- hist(x, breaks = breaks, plot = FALSE)
  data.frame(group = group, xmin = head(h$breaks, -1), xmax = h$breaks[-1],
             prop = h$counts / length(x))
}

# Shared look of the get_fig_eval() panels
theme_fig_eval <- function() {
  theme_minimal(base_size = 11) +
    theme(
      text = element_text(colour = "#0b0b0b"),
      axis.text = element_text(colour = "#52514e"),
      axis.title = element_text(colour = "#52514e"),
      plot.title = element_text(face = "bold", size = rel(1)),
      plot.title.position = "plot",
      panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(colour = "#52514e"),
      legend.title = element_text(colour = "#52514e"),
      legend.key.size = grid::unit(0.9, "lines"),
      legend.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 14, 8, 8)
    )
}

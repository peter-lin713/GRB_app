#' mc_plots.R — diagnostic plots for Monte Carlo redshift predictions.
#'
#' Turns the MC summary produced by summarize_mc() (see mc_error_propagation.R)
#' into predicted-vs-observed scatter plots with 95% MC error bars. Points are
#' coloured by whether the truth falls inside each GRB's MC interval ("inside
#' the cone"), which is the headline calibration check.
#'
#' Helpers (top to bottom):
#'   inside_mc_cone     logical flag: does the truth lie in the MC interval?
#'   pred_vs_obs_plot   shared ggplot builder (errorbars + identity line)
#'   metrics_subtitle   "n / r / RMSE" string for plot subtitles
#'   make_mc_plots      writes the three diagnostic PNGs (entry point)
#'
#' make_mc_plots() is the public entry point called by the pipeline.

library(ggplot2)

#' Flag whether each GRB's MC 95% interval covers the truth ("inside cone").
#'
#' @param mc_summary Dataframe from summarize_mc() (has log10z_lower/upper)
#' @param y_true Named numeric vector of truths, names match rownames(mc_summary)
#' @return Logical vector of length nrow(mc_summary)
inside_mc_cone <- function(mc_summary, y_true) {
  y <- y_true[rownames(mc_summary)]
  y >= mc_summary$log10z_lower & y <= mc_summary$log10z_upper
}

#' Common predicted-vs-observed plot with MC error bars.
#' @param df Dataframe with columns x, y, y_lo, y_hi, inside (logical)
#' @param title Plot title
#' @param xlab,ylab Axis labels
#' @param subtitle Optional subtitle (e.g. metrics string)
#' @param show_outliers If FALSE, drop points where inside == FALSE
pred_vs_obs_plot <- function(df, title, xlab, ylab, subtitle = NULL,
                             show_outliers = TRUE) {
  if (!show_outliers) df <- df[df$inside, ]

  lims <- range(c(df$x, df$y, df$y_lo, df$y_hi), finite = TRUE)

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray40") +
    geom_errorbar(aes(ymin = y_lo, ymax = y_hi, color = inside),
                  width = 0, alpha = 0.5) +
    geom_point(aes(color = inside, shape = inside), size = 1.8) +
    scale_color_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "#d7301f"),
                       labels = c(`TRUE` = "Inside 95% MC cone",
                                  `FALSE` = "Outside 95% MC cone"),
                       name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                       labels = c(`TRUE` = "Inside 95% MC cone",
                                  `FALSE` = "Outside 95% MC cone"),
                       name = NULL) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
  p
}

#' Metrics string for plot subtitle.
metrics_subtitle <- function(x, y) {
  n   <- length(x)
  r   <- cor(x, y)
  rms <- sqrt(mean((x - y)^2))
  sprintf("n = %d   |   r = %.3f   |   RMSE = %.3f", n, r, rms)
}

#' Generate and save the three MC diagnostic plots.
#'
#' @param mc_summary Dataframe from summarize_mc(), indexed by GRB
#' @param y_true Named vector of true log10(z+1) values
#' @param z_true Named vector of true z values (linear)
#' @param out_dir Directory to write PNGs into (trailing slash ok)
#' @return (invisibly) a list of the three ggplot objects
make_mc_plots <- function(mc_summary, y_true, z_true, out_dir) {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Align truths to summary rows
  grbs <- rownames(mc_summary)
  yt   <- y_true[grbs]
  zt   <- z_true[grbs]
  inside <- inside_mc_cone(mc_summary, y_true)

  # ----- Plot 1: all points, log10 scale -----
  df1 <- data.frame(
    x      = yt,
    y      = mc_summary$log10z_mean,
    y_lo   = mc_summary$log10z_lower,
    y_hi   = mc_summary$log10z_upper,
    inside = inside
  )
  sub1 <- sprintf("%s   |   %.1f%% inside cone",
                  metrics_subtitle(df1$x, df1$y),
                  100 * mean(inside))
  p1 <- pred_vs_obs_plot(
    df1,
    title = "Predicted vs observed log10(z+1) — all points",
    subtitle = sub1,
    xlab = "Observed log10(z+1)",
    ylab = "Predicted log10(z+1)  (MC mean, 95% interval)"
  )
  ggsave(file.path(out_dir, "mc_pred_vs_obs_log10_all.png"),
         p1, width = 6.5, height = 6.5, dpi = 300)

  # ----- Plot 2: inside-cone only, log10 scale -----
  df2 <- df1[df1$inside, ]
  sub2 <- metrics_subtitle(df2$x, df2$y)
  p2 <- pred_vs_obs_plot(
    df2,
    title = "Predicted vs observed log10(z+1) — inside 95% MC cone",
    subtitle = sub2,
    xlab = "Observed log10(z+1)",
    ylab = "Predicted log10(z+1)  (MC mean, 95% interval)",
    show_outliers = FALSE
  )
  ggsave(file.path(out_dir, "mc_pred_vs_obs_log10_inside.png"),
         p2, width = 6.5, height = 6.5, dpi = 300)

  # ----- Plot 3: linear z, inside-cone only -----
  df3 <- data.frame(
    x      = zt[inside],
    y      = mc_summary$z_mean[inside],
    y_lo   = mc_summary$z_lower[inside],
    y_hi   = mc_summary$z_upper[inside],
    inside = TRUE
  )
  sub3 <- metrics_subtitle(df3$x, df3$y)
  p3 <- pred_vs_obs_plot(
    df3,
    title = "Predicted vs observed z (linear) — inside 95% MC cone",
    subtitle = sub3,
    xlab = "Observed z",
    ylab = "Predicted z  (MC mean, 95% interval)",
    show_outliers = FALSE
  )
  ggsave(file.path(out_dir, "mc_pred_vs_obs_z_inside.png"),
         p3, width = 6.5, height = 6.5, dpi = 300)

  invisible(list(all_log10 = p1, inside_log10 = p2, inside_z = p3))
}

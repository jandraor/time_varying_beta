plot_ll_se <- function(sens_df, clr) {
  
  ggplot(sens_df, aes(x = Np / 1000, y = loglik.se)) + 
    geom_point(colour = clr) +
    geom_line(alpha = 0.15, colour = clr, size = 1) +
    geom_text_repel(aes(label = round(loglik, 0)), size = 2, colour = "grey55") +
    scale_y_log10() +
    scale_x_continuous(labels = comma) +
    geom_hline(yintercept = 2, colour = "grey85", linetype = "dotted") +
    facet_grid(inv_dt~id) +
    labs(x = "Number of particles (Thousands)", y = "Log-lik SE (Log scale)",
         caption  = "Labels within plots indicate log-lik values",
         title    = "Sensitivity of the Log-lik standard error (SE)",
         subtitle = "Testing points (columns) vs intervals per dt (rows) ") +
    theme_bw() +
    theme(panel.grid    = element_blank(),
          axis.title = element_text(colour = "grey40"),
          axis.text     = element_text(colour = "grey60", size = 6),
          plot.title    = element_text(colour = "grey25", size = 10),
          plot.subtitle = element_text(colour = "grey25", size = 8),
          axis.ticks    = element_line(colour = "grey60", size = 0.25),
          axis.line     = element_line(colour = "grey60", size = 0.25))
}

plot_ll_time <- function(sens_df, clr) {
  
  ggplot(sens_df, aes(x = Np / 1000, y = time)) + 
    geom_point(colour = clr, aes(shape = as.factor(inv_dt))) +
    geom_line(aes(group = inv_dt),alpha = 0.15, colour = clr, size = 1) +
    scale_y_log10() +
    scale_x_continuous(labels = comma) +
    labs(x = "Number of particles (Thousands)", y = "Minutes (Log scale)",
         title    = "Time sensitivity on testing point 1",
         shape = "# of intervals")  +
    theme_pubr() +
    theme(axis.title    = element_text(colour = "grey40"),
          axis.text     = element_text(colour = "grey60", size = 6),
          plot.title    = element_text(colour = "grey25", size = 10),
          axis.ticks    = element_line(colour = "grey60", size = 0.25),
          axis.line     = element_line(colour = "grey60", size = 0.25))
}


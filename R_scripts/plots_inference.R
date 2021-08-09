par_comparison_stochastic <- function(df) {
  ggplot(df, aes(x = method, y = mle)) +
    geom_errorbar(aes(ymin = ll, ymax = ul, colour = method)) +
    geom_point(aes(colour = method)) +
    scale_colour_manual(values = c("#1D2F31", "#136AA0", "#C1A284")) +
    facet_wrap(~name, labeller = label_parsed, scales = "free_y") +
    theme_pubr() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 6))
}

plot_wkl_fit <- function(summary_fit, actual_df, y_lab, title_lab, shape, 
                         plot_colour) {
  
  ggplot(summary_fit, aes(x = week, y = q50)) +
    geom_line(colour = plot_colour) +
    scale_y_continuous(labels = comma) +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5),
                fill = plot_colour) +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75),
                fill = plot_colour) +
    geom_point(data = actual_df, aes(x = week, y = y), colour = data_colour,
               size = 2, shape = shape, alpha = 0.95) +
    scale_x_continuous(breaks = 1:11) +
    theme_pubr() +
    labs(y = parse(text = y_lab), x = "Week", title = title_lab) +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25))
  
}

plot_hidden_re <- function(re_df, pred_colour) {
  ggplot(re_df, aes(x = week, y = q50)) +
    geom_line(colour = pred_colour) +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5),
                fill = pred_colour) +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75),
                fill = pred_colour) +
    geom_vline(xintercept = 1.86, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 1.96, y = 6, size = 1.5,
             label = "Delay phase", hjust = 0, size = 3, colour = "grey50") +
    geom_vline(xintercept = 4, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 4.1, y = 6, size = 1.5,
             label = "Stay at home", hjust = 0, size = 3, colour = "grey50") +
    geom_hline(yintercept = 1, colour = ET_colour, linetype = "dashed") +
    labs(y = parse(text = "Re[t]"), x = "Week",
         title = "C) Effective reproductive number")+
    scale_x_continuous(breaks = 1:11) +
    theme_pubr()  +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25))
}

plot_filt_dist_re <- function(samples_df, wkl_inc, dist_col) {
  
  ggplot(samples_df, aes(x = as.factor(week), y = value)) +
    geom_violin(colour = dist_col) +
    facet_wrap(~ as.factor(week), scales = "free") +
    geom_hline(data = wkl_inc, aes(x = as.factor(week)), yintercept = 1,
               linetype = "dashed", colour = ET_colour) +
    theme_classic() +
    labs(x = "Week", y = parse(text = "Re[t]")) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y  = element_text(size = 4),
          axis.line  = element_line(colour = "grey60", size = 0.2),
          axis.ticks = element_line(colour = "grey60", size = 0.2),
          axis.title = element_text(size = 8, colour = "grey40"))
}

plot_filt_dist <- function(samples_df, wkl_data, y_label, dist_col) {
  
  ggplot(samples_df, aes(x = as.factor(week), y = value, group = week)) +
    geom_violin(colour = dist_col) +
    facet_wrap(~time, scales = "free") +
    geom_hline(data = wkl_data, aes(x = as.factor(week), yintercept = y), 
               linetype = "dotted", colour = data_colour) +
    theme_classic() +
    labs(x = "Week", y = parse(text = y_label)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y  = element_text(size = 4),
          axis.line  = element_line(colour = "grey60", size = 0.2),
          axis.ticks = element_line(colour = "grey60", size = 0.2),
          axis.title = element_text(size = 8, colour = "grey40"))
}

plot_fits_by_order <- function(sim_df, actual_data, df_labels, shape, x_lab,
                               y_lab) {
  ggplot(sim_df, aes(x = time, y = value)) +
    geom_line(aes(group = iter), alpha = 0.25, colour = STH_colour) +
    geom_point(data = actual_data, aes(y = y), colour = data_colour,
               alpha = 0.8, size = 0.25, shape = shape) +
    facet_wrap(~order, ncol = 1) +
    geom_text(data = df_labels, aes(label = label, x = x, y = y), 
                       size = 2, colour = "grey50") +
    labs(y = parse(text = y_lab), x = x_lab) +
    theme_test() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 6, colour = "grey60"),
          axis.line  = element_line(colour = "grey60", size = 0.25),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.title = element_text(size = 8, colour = "grey40"),
          panel.border = element_rect(colour = "grey95", fill=NA, size = 0.5))
}

plot_fit_comparison <- function(sim_data, actual_data, y_label, title_label, 
                                shape){
  ggplot(sim_data, aes(x = week, y = q50)) +
    geom_line(aes(colour = DGP, group = DGP)) +
    scale_colour_manual(values = sim_colours) +
    scale_fill_manual(values = sim_colours) +
    scale_y_continuous(labels = comma) +
    geom_ribbon(alpha = 0.25, aes(ymin = q2.5, ymax = q97.5,
                                  fill = DGP)) +
    geom_point(data = actual_data, aes(x = week, y = y), colour = data_colour,
               size = 1, shape = shape) +
    theme_pubr() +
    labs(y = parse(text = y_label), x = "Week", 
         title = title_label) +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(size = 9, colour = "grey25"),
          axis.ticks = element_line(colour = "grey60"),
          axis.line  = element_line(colour = "grey60"))
}




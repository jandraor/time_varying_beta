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

plot_wkl_fit <- function(summary_fit, actual_df, y_lab, title_lab, shape) {
  
  ggplot(summary_fit, aes(x = week, y = q50)) +
    geom_line(colour = "#C7C2F9") +
    scale_y_continuous(labels = comma) +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5),
                fill = "#C7C2F9") +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75),
                fill = "#C7C2F9") +
    geom_point(data = actual_df, aes(x = week, y = y), colour = "#194973",
               size = 2, shape = shape, alpha = 0.95) +
    scale_x_continuous(breaks = 1:11) +
    theme_pubr() +
    labs(y = y_lab, x = "Week", title = title_lab) +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(size = 9, colour = "grey25"),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25))
  
}

plot_hidden_re <- function(re_df) {
  ggplot(re_df, aes(x = week, y = q50)) +
    geom_line(colour = "#C7C2F9") +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5),
                fill = "#C7C2F9") +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75),
                fill = "#C7C2F9") +
    geom_vline(xintercept = 1.86, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 1.96, y = 6, size = 1.5,
             label = "Delay phase", hjust = 0, size = 3, colour = "grey50") +
    geom_vline(xintercept = 4, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 4.1, y = 6, size = 1.5,
             label = "Stay at home", hjust = 0, size = 3, colour = "grey50") +
    geom_hline(yintercept = 1, colour = "#F9773B", linetype = "dashed") +
    labs(y = parse(text = "R[e]"), x = "Week",
         title = "C) Predicted effective reproductive number")+
    scale_x_continuous(breaks = 1:11) +
    theme_pubr()  +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(size = 9, colour = "grey25"),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25))
}

plot_filt_dist_re <- function(samples_df, wkl_inc) {
  
  ggplot(samples_df, aes(x = as.factor(week), y = value)) +
    geom_violin(colour = "#C7C2F9") +
    facet_wrap(~ as.factor(week), scales = "free") +
    geom_hline(data = wkl_inc, aes(x = as.factor(week)), yintercept = 1,
               linetype = "dashed", colour = "#F9773B") +
    theme_classic() +
    labs(x = "Week", y = parse(text = "R[e]")) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y  = element_text(size = 4),
          axis.line  = element_line(colour = "grey60", size = 0.2),
          axis.ticks = element_line(colour = "grey60", size = 0.2),
          axis.title = element_text(size = 8, colour = "grey40"))
}

plot_filt_dist <- function(samples_df, wkl_data, y_label) {
  
  ggplot(samples_df, aes(x = as.factor(week), y = value, group = week)) +
    geom_violin(colour = "#C7C2F9") +
    facet_wrap(~time, scales = "free") +
    geom_hline(data = wkl_data, aes(x = as.factor(week), yintercept = y), 
               linetype = "dotted", colour = "#194973") +
    theme_classic() +
    labs(x = "Week", y = y_label) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y  = element_text(size = 4),
          axis.line  = element_line(colour = "grey60", size = 0.2),
          axis.ticks = element_line(colour = "grey60", size = 0.2),
          axis.title = element_text(size = 8, colour = "grey40"))
}



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

plot_hidden_re <- function(re_df) {
  ggplot(re_df, aes(x = week, y = median)) +
    geom_hline(yintercept = 1, colour = "grey60", linetype = "dashed") +
    geom_line(colour = "steelblue") +
    geom_ribbon(alpha = 0.25, aes(ymin = lower_lim, ymax = upper_lim),
                fill = "steelblue") +
    geom_vline(xintercept = 1.86, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 1.96, y = 6, size = 1.5,
             label = "Delay phase", hjust = 0, size = 3, colour = "grey50") +
    geom_vline(xintercept = 4, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 4.1, y = 6, size = 1.5,
             label = "Stay at home", hjust = 0, size = 3, colour = "grey50") +
    labs(y = parse(text = "Re[t]"), x = "Week",
         title = "Effective reproductive number")+
    scale_x_continuous(breaks = 1:11) +
    theme_pubr()  +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(size = 9, colour = "grey25"),
          axis.ticks = element_line(colour = "grey60"))
}

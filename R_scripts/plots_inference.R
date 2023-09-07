par_comparison_stochastic <- function(df) {
  
  df <- df |>  
    mutate(name = case_when(name == "P_0" ~ "P[0]",
                            name == "Re_0" ~ '\u211c[0]',
                            TRUE ~ as.character(name)))
  
  ggplot(df, aes(x = method, y = mle)) +
    geom_errorbar(aes(ymin = ll, ymax = ul, colour = method)) +
    geom_point(aes(colour = method)) +
    scale_colour_manual(values = c("#1D2F31", "#136AA0", "#C1A284")) +
    facet_wrap(~name, labeller = label_parsed, scales = "free_y") +
    labs(x = "Method", y = "Value") +
    theme_pubr() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 6),
          text = element_text(family = "Arial Unicode MS"))
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
    scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
    theme_pubr() +
    labs(y = parse(text = y_lab), x = "Week", title = title_lab) +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25),
          text = element_text(family = "Arial Unicode MS"))
  
}

plot_hidden_re <- function(re_df, pred_colour, plot_title) {
  
  ggplot(re_df, aes(x = week, y = q50)) +
    geom_line(colour = pred_colour) +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5),
                fill = pred_colour) +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75),
                fill = pred_colour) +
    geom_vline(xintercept = 1.86, colour = "grey50", linetype = "dotted", 
               alpha = 0.75) +
    annotate("text", x = 1.96, y = 6, size = 1.25,
             label = "Delay phase", hjust = 0, size = 1.5, colour = "grey50") +
    geom_vline(xintercept = 4, colour = "grey50", linetype = "dotted",
               alpha = 0.75) +
    annotate("text", x = 4.1, y = 6.5, size = 1.25,
             label = "Stay at home", hjust = 0, size = 1, colour = "grey50") +
    geom_hline(yintercept = 1, colour = ET_colour, linetype = "dashed") +
    labs(y = parse(text = "\u211c[t]"), x = "Week",
         title = plot_title)+
    scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
    theme_pubr()  +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25),
          text = element_text(family = "Arial Unicode MS"))
}

plot_filt_dist_re <- function(samples_df, wkl_inc, dist_col) {
  
  ggplot(samples_df, aes(x = as.factor(week), y = value)) +
    geom_violin(colour = dist_col) +
    facet_wrap(~ as.factor(week), scales = "free") +
    geom_hline(data = wkl_inc, aes(x = as.factor(week)), yintercept = 1,
               linetype = "dashed", colour = ET_colour) +
    scale_y_continuous(n.breaks = 4) +
    theme_classic() +
    labs(x = "Week", y = parse(text = "\u211c[t]")) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y  = element_text(size = 5),
          axis.line  = element_line(colour = "grey60", size = 0.2),
          axis.ticks = element_line(colour = "grey60", size = 0.2),
          axis.title = element_text(size = 8, colour = "grey40"),
          text = element_text(family = "Arial Unicode MS"))
}

plot_daily_re <- function(re_df, pred_colour, plot_title) {
  
  ggplot(re_df, aes(x = time, y = q50)) +
    geom_line(colour = pred_colour) +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5),
                fill = pred_colour) +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75),
                fill = pred_colour) +
    geom_hline(yintercept = 1, colour = ET_colour, linetype = "dashed") +
    labs(y = parse(text = "\u211c[t]"), x = "Day",
         title = plot_title)+
    theme_pubr()  +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25),
          text = element_text(family = "Arial Unicode MS"))
}

plot_filt_dist <- function(samples_df, wkl_data, y_label, dist_col) {
  
  ggplot(samples_df, aes(x = as.factor(week), y = value)) +
    geom_violin(colour = dist_col) +
    facet_wrap(~as.factor(week), scales = "free") +
    geom_hline(data = wkl_data, aes(x = as.factor(week), yintercept = y), 
               linetype = "dotted", colour = data_colour) +
    scale_y_continuous(n.breaks = 4) +
    theme_classic() +
    labs(x = "Week", y = parse(text = y_label)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y  = element_text(size = 5),
          axis.line  = element_line(colour = "grey60", size = 0.2),
          axis.ticks = element_line(colour = "grey60", size = 0.2),
          axis.title = element_text(size = 8, colour = "grey40"),
          text = element_text(family = "Arial Unicode MS"))
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
          panel.border = element_rect(colour = "grey95", fill = NA, size = 0.5))
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
          axis.line  = element_line(colour = "grey60"),
          text = element_text(family = "Arial Unicode MS"))
}

plot_MLE_neighbourhood <- function(df, col) {
  
  df <- rename(df, `P[0]` = P_0)
  
  ggpairs(df, 
          labeller = label_parsed, 
          diag = list(continuous = wrap('densityDiag', colour = col)),
          lower = list(continuous = wrap("points", colour = col))) +
    theme_pubclean() +
    theme(axis.text = element_text(size = 6))

}

plot_daily_fit <- function(pred_df, data_df, y_lab, title_lab, shape, 
                           plot_colour, point_size = 2) {
  
  ggplot(pred_df, aes(x = time, y = q50)) +
    geom_line(colour = plot_colour) +
    scale_y_continuous(labels = comma) +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5),
                fill = plot_colour) +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75),
                fill = plot_colour) +
    geom_point(data = data_df, aes(y = y), colour = data_colour,
               size = point_size, shape = shape, alpha = 0.95) +
    theme_pubr() +
    labs(x = "Day", y = parse(text = y_lab), title = title_lab) +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25))
}

plot_fit_by_chain_type <- function(df, irish_data, order) {
  
  n_iters      <- unique(df$iter)
  sample_iters <- sample(n_iters, 100, replace = FALSE)
  df           <- df %>% filter(iter %in% sample_iters)
  
  ggplot(df, aes(x = time, y = value)) +
    geom_line(aes(group = iter, colour = fit), alpha = 1) +
    geom_point(data = irish_data, aes(y = y), colour = data_colour,
               alpha = 0.8, size = 0.5, shape = 16) +
    facet_wrap(~fit, scales = "free") +
    scale_colour_manual(values = c("#007C7C", "#bcdcdc")) +
    theme_pubr() +
    labs(title    = str_glue("Order: {order}"),
         subtitle = str_glue("Do the samples fit the incidence data?"),
         x = "", 
         y = parse(text = "C[t]")) +
    theme(axis.text = element_text(colour = "grey60", size = 8),
          axis.line  = element_line(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          axis.title = element_text(colour = "grey40"),
          legend.position = "none")
}  

plot_ll_by_order <- function(lik_df) {
  ggplot(lik_df, aes(x = as.factor(order), y = log_lik)) +
    geom_violin(colour = STH_colour) +
    theme_pubr() +
    labs(x = "Order", y = "Log-likelihood")
}

plot_bp_by_par <- function(summary_var, var_name) {
  
  var_df <- summary_var |> filter(var == var_name)
  
  ggplot(var_df, aes(x = as.factor(chain), y = value)) +
    geom_boxplot(aes(colour = converges)) +
    scale_colour_manual(values = c("#007C7C", "#bcdcdc")) +
    theme_classic() +
    facet_wrap(~order, scales = "free") +
    labs(x = "Chain", y = "Value", 
         subtitle = parse(text = str_glue('"Variable: " ~{var_name}'))) +
    theme(axis.text = element_text(colour = "grey60", size = 8),
          axis.line  = element_line(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          legend.position = "top")
}

plot_sth_wkl_rt <- function(summary_Re) {
  
  ggplot(summary_Re, aes(x = week, y = q50)) +
    geom_hline(yintercept = 1, colour = ET_colour, linetype = "dashed") +
    geom_line(colour = STH_colour) +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5),
                fill = STH_colour) +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75),
                fill = STH_colour) +
    geom_vline(xintercept = 1.86, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 1.96, y = 6, size = 1.5,
             label = "Delay phase", hjust = 0, size = 3, colour = "grey50") +
    geom_vline(xintercept = 4, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 4.1, y = 6, size = 1.5,
             label = "Stay at home", hjust = 0, size = 3, colour = "grey50") +
    labs(y = parse(text = "Re[t]"), x = "Week",
         title = "Effective reproductive number")+
    scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
    scale_y_continuous(limits = c(0, 8)) +
    theme_pubr()  +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25))
}

plot_par_comparison <- function(par_summary) {
  
  par_summary <- par_summary |> 
    mutate(name = case_when(name == "P_0" ~ "P[0]",
                            name == "Re_0" ~ '\u211c[0]',
                            TRUE ~ as.character(name)))
  
  ggplot(par_summary, aes(x = model, y = mle)) +
    geom_errorbar(aes(ymin = ll, ymax = ul), colour = GBM_colour) +
    geom_point(colour = GBM_colour) +
    facet_wrap(~name, labeller = label_parsed, scales = "free") +
    theme_pubr() +
    labs(x = "") +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25),
          text = element_text(family = "Arial Unicode MS"))
}

plot_daily_fit_by_order <- function(pred_list, data_df, data_shape, y_label) {
  
  plot_list <- imap(pred_list, function(df, i) {
    
    iters <- unique(df$iter)
    
    sample_iters <- sample(iters, 100, replace = FALSE)
    df           <- df |>  filter(iter %in% sample_iters)
    
    ggplot(df, aes(x = time, y = value)) +
      geom_line(aes(group = iter), colour = STH_colour, alpha = 0.5) +
      theme_pubr() +
      geom_point(data = data_df, aes(y = y), colour = data_colour, 
                 alpha = 0.8, size = 0.5, shape = data_shape) +
      labs(subtitle = str_glue("Order: {i}"),
           x = "", y = parse(text = y_label)) +
      theme(axis.text = element_text(colour = "grey60", size = 8),
            axis.line  = element_line(colour = "grey60"),
            axis.ticks = element_line(colour = "grey60"),
            axis.title = element_text(colour = "grey40")) 
  }) 
  
  wrap_plots(plot_list) +
    plot_annotation(caption = "Time since the first reported case [Days]",
                    theme = theme(plot.caption = element_text(hjust = 0.5)))
}

plot_daily_fit_comparison <- function(pois_df, nbin_df, data_df, order, 
                                      data_shape, y_label) {
  
  iters_pois   <- unique(pois_df$iter)
  sample_iters <- sample(iters_pois, 100, replace = FALSE)
  pois_df      <- pois_df |>  filter(iter %in% sample_iters) |> 
    mutate(dist = "Pois")
  
  
  iters_nbin   <- unique(nbin_df$iter)
  sample_iters <- sample(iters_nbin, 100, replace = FALSE)
  nbin_df      <- nbin_df |>  filter(iter %in% sample_iters) |> 
    mutate(dist = "Nbin")
  
  df <- bind_rows(pois_df, nbin_df)
  
  ggplot(df, aes(x = time, y = value)) +
    geom_line(aes(group = iter), colour = STH_colour, alpha = 0.5) +
    facet_wrap(~dist) +
    theme_pubr() +
    geom_point(data = data_df, aes(y = y), colour = data_colour, 
               alpha = 0.8, size = 0.5, shape = data_shape) +
    labs(subtitle = str_glue("Order: {order}"),
         x = "", y = parse(text = y_label)) +
    theme(axis.text = element_text(colour = "grey60", size = 8),
          axis.line  = element_line(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          axis.title = element_text(colour = "grey40")) 
}

plot_re_comparison <- function(df) {
  ggplot(df, aes(x = week, y = q50)) +
    geom_hline(yintercept = 1, colour = ET_colour, linetype = "dashed") +
    geom_line(aes(group = dist, colour = dist)) +
    geom_ribbon(alpha = 0.2, aes(ymin = q2.5, ymax = q97.5, group = dist,
                                 fill = dist)) +
    geom_ribbon(alpha = 0.5, aes(ymin = q25, ymax = q75, group = dist, 
                                 fill = dist)) +
    scale_fill_manual(values = c(STH_colour, "grey50")) +
    scale_colour_manual(values = c(STH_colour, "grey50")) +
    geom_vline(xintercept = 1.86, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 1.96, y = 6, size = 1.5,
             label = "Delay phase", hjust = 0, size = 3, colour = "grey50") +
    geom_vline(xintercept = 4, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 4.1, y = 6, size = 1.5,
             label = "Stay at home", hjust = 0, size = 3, colour = "grey50") +
    labs(y = parse(text = "\u211c[t]"), x = "Week",
         title = "Effective reproductive number")+
    scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
    scale_y_continuous(limits = c(0, 8)) +
    theme_pubr()  +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(colour = "grey25", size = 8),
          axis.ticks = element_line(colour = "grey60", size = 0.25),
          axis.line  = element_line(colour = "grey60", size = 0.25),
          text = element_text(family = "Arial Unicode MS"))
}

plot_re_by_DGP <- function (df) {
  
  ggplot(df, aes(x = week, y = q50)) +
    geom_hline(yintercept = 1, colour = ET_colour, linetype = "dashed") +
    geom_line(aes(group = DGP, colour = DGP)) +
    geom_ribbon(alpha = 0.25, aes(ymin = q2.5, ymax = q97.5,
                                  fill = DGP)) +
    geom_vline(xintercept = 1.86, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 1.96, y = 6, size = 1.5,
             label = "Delay phase", hjust = 0, size = 3, colour = "grey50") +
    geom_vline(xintercept = 4, colour = "grey50", linetype = "dotted") +
    annotate("text", x = 4.1, y = 6, size = 1.5,
             label = "Stay at home", hjust = 0, size = 3, colour = "grey50") +
    scale_colour_manual(values = sim_colours) +
    scale_fill_manual(values = sim_colours) +
    labs(y = parse(text = "\u211c[t]"), x = "Week",
         title = "C) Effective reproductive number per DGP")+
    scale_x_continuous(breaks = 1:11) +
    theme_pubr()  +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(size = 9, colour = "grey25"),
          axis.ticks = element_line(colour = "grey60"),
          axis.line  = element_line(colour = "grey60"),
          text = element_text(family = "Arial Unicode MS"))
}

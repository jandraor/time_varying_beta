library(extraDistr)
library(GGally)
library(ggalt)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(scales)

source("./R_scripts/plots_inference.R")

data_colour <- "#FF7241"
GBM_colour  <- "#344D77"
CIR_colour  <- "#C55A82"
STH_colour  <- "#62DCDA"
ET_colour   <- "#FFD34B"

sim_colours <- c(GBM_colour, CIR_colour, STH_colour)

pm_likelihood <- function(df) {
  tidy_loglik_df <- df %>% filter(type != "guess") %>%
    mutate(row_number = row_number()) %>% 
    pivot_longer(c(tau, sigma), names_to = "var", values_to = "var_value")
  
  tidy_loglik_df %>% filter(loglik > max(loglik)-10) %>%
    ggplot(aes(x = var_value, y = loglik))+
    geom_point(colour = "grey40")+
    facet_wrap(~var, scales = "free_x", labeller = label_parsed) +
    theme_classic2() +
    labs(
      x     = "Value",
      title ="Profile likelihood"
    )
}

hist_final_epi_size <- function(sim_tc_df) {
  
  ggplot(sim_tc_df, aes(x = magnitude_order)) +
    geom_histogram(stat = "count", aes(fill = colour)) +
    scale_fill_manual(values = c("grey60", "steelblue")) +
    labs(y = "Count", x = "Final epidemic size") +
    theme_pubclean() +
    theme(legend.position = "none")
}

daily_epicurve <- function(irish_data, formatted_tc) {
  ggplot(irish_data, aes(x = date, y = y)) +
    geom_col(colour = "white", fill = data_colour, alpha = 0.95) +
    theme_pubr() +
    geom_vline(xintercept = as_date("2020-02-29") + 13, colour = "grey50",
               linetype = "dotted") +
    geom_vline(xintercept = as_date("2020-02-29") + 28, colour = "grey10",
               linetype = "dotted") +
    annotate("text", x = as_date("2020-02-29") + 13.5, y = 750, size = 1.5,
             label = "Delay phase", hjust = 0, size = 3, colour = "grey50") +
    annotate("text", x = as_date("2020-02-29") + 28.5, y = 920, size = 1.5,
             label = "Stay at home", hjust = 0, size = 3, colour = "grey50") +
    labs(x     = "Date of Lab specimen collection",
         y     = "Incidence [new cases per day]",
         title = "Daily detected cases",
         subtitle = str_glue("Total cases: {formatted_tc}")) +
    theme(plot.subtitle = element_text(size = 8, colour = "grey60"),
          plot.title = element_text(size = 9, colour = "grey10"),
          axis.title = element_text(size = 8, colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60", size = 6),
          axis.ticks = element_line(colour = "grey60"))
}

daily_epi_trend <- function(irish_data, title) {
  
  ggplot(irish_data, aes(x = date, y = y)) +
    geom_vline(xintercept = as_date("2020-02-29") + 13, colour = "grey50",
               linetype = "dotted") +
    geom_vline(xintercept = as_date("2020-02-29") + 28, colour = "grey10",
               linetype = "dotted") +
    annotate("text", x = as_date("2020-02-29") + 13.5, y = 750, size = 1.5,
             label = "Delay phase", hjust = 0, size = 3, colour = "grey50") +
    annotate("text", x = as_date("2020-02-29") + 28.5, y = 930, size = 1.5,
             label = "Stay at home", hjust = 0, size = 3, colour = "grey50") +
    labs(x = "Date of Lab specimen collection",
         y = "Incidence [new cases per day]",
         title = title) +
    geom_point(colour = data_colour, size = 0.75, shape = 18, alpha = 0.95) +
    geom_line(stat="smooth", method = "loess", formula = y ~ x, alpha = 0.25,
              colour = data_colour, span = 0.25, size = 1) +
    theme_pubr() +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(size = 9, colour = "grey10"),
          axis.ticks = element_line(colour = "grey60"),
          axis.line =  element_line(colour = "grey80"))
}

weekly_epicurve <- function(wkl_data, title) {
  ggplot(wkl_data, aes(x = week, y = y)) +
    geom_col(fill = data_colour) +
    scale_x_continuous(breaks = 0:12) +
    theme_pubclean() +
    labs(y = "Incidence [new cases per week]",
         x = "Week of Lab specimen collection",
         title = title) +
  theme_pubr() +
    theme(plot.title = element_text(size = 9),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60", size = 7),
          axis.ticks = element_line(colour = "grey60"),
          axis.title = element_text(size = 8, colour = "grey40"))
}

incidence_graphs <- function(irish_data, formatted_tc, wkl_inc) {
  g1 <- daily_epicurve(irish_data, formatted_tc)
  g2 <- daily_epi_trend(irish_data)
  g3 <- weekly_epicurve(wkl_inc)
  
  (g1/g2) | g3
}

plot_imputed_mob <- function(data, raw_data, imputed_data) {
  
  imp_df <- data.frame(time         = dates,
                       raw_data     = raw_data,
                       imputed_data = imputed_data) %>% 
    mutate(is.imp = ifelse(is.na(raw_data), "imputed values", "known values"))
  
  ggplot(imp_df, aes(x = time, y = imputed_data)) +
    geom_line(colour = data_colour, alpha = 0.3) +
    geom_point(shape = 16, aes(colour = is.imp)) +
    scale_colour_manual(values = c("blue", data_colour), name = "") +
    theme_pubr() +
    labs(x = "Date", y = "Index",
         title = "Ireland's driving data",
         subtitle = "Missing data has been replaced by linear interpolation estimates",
         caption = "This dataset has been normalised by its initial value (141)") +
    theme(plot.subtitle = element_text(colour = "grey65"),
          plot.caption  = element_text(colour = "grey65"),
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"))
}

plot_daily_mobility <- function(daily_mob_df) {
  
  daily_mob_df <- daily_mob_df %>% 
    mutate(time        = row_number(),
           end_of_week = ifelse(time %% 7 == 0, TRUE, FALSE)) %>% 
    filter(time <= 77)
  
  ggplot(daily_mob_df, aes(x = date, y = y)) +
    geom_point(colour = data_colour, aes(alpha = end_of_week), size = 0.5)   +
    scale_alpha_manual(values = c(0.25, 1)) +
    geom_line(stat = "smooth", method = "loess", formula = y ~ x, alpha = 0.1,
              colour = data_colour, span = 0.25, size = 1) +
    labs(x = "Date", y = "Mobility index",
         title = "B) Mobility trend (Driving)") +
    theme_pubr() +
    theme(legend.position = "none",
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(size = 9, colour = "grey10"),
          axis.ticks = element_line(colour = "grey60"),
          axis.line =  element_line(colour = "grey60"),
          axis.title = element_text(size = 8, colour = "grey40"))
  
}

plot_daily_data <- function(irish_data, drv_df) {
  
  g1 <- daily_epi_trend(irish_data)
  
  g2 <- plot_daily_mobility(daily_mob_df) +
    geom_point(colour = data_colour, size = 0.5)
  
  g1 / g2
}

loglik_traces <- function(traces_df, lims, traces_col) {
  
  traces_df %>% filter(variable == "loglik") %>% 
    ggplot(aes(x = iteration,y = value,group = L1))+
    geom_line(colour = traces_col, alpha = 0.25)+
    scale_y_continuous(limits = lims) +
    guides(color=FALSE) +
    labs(x = "Iteration", y = "Log-likelihood") +
    theme_pubclean()

  }

loglik_se_sens <- function(sens_results, DGP_colour) {
  ggplot(sens_results, aes(x = Np, y = loglik.se)) +
    scale_x_log10(labels = comma) +
    geom_line(stat = "smooth", method = "loess", formula = y ~ x, alpha = 0.5,
              colour = DGP_colour, span = 0.75, size = 1, linetype = "dashed") + 
    geom_point(colour = DGP_colour, size = 4) +
    theme_pubclean() +
    labs(y = "Log-likelihood SE", x = "Number of particles (Log scale)")
}

mif_traces <- function(traces_df, unknown_pars, traces_col) {
  traces_df %>% 
    filter(variable %in% c("loglik", unknown_pars)) %>%
    ggplot(aes(x= iteration,y = value, group = L1))+
    geom_line(color = traces_col, alpha = 0.25) +
    theme_pubr() +
    scale_y_continuous(labels = comma) +
    labs(y = "Value", x = "Iteration") +
    facet_wrap(~variable,scales = "free_y", labeller = label_parsed)
}

raw_likelihood <- function(all_df, var_x, point_colour) {
  
  var_x <- as.character(var_x)
  lab_x <- var_x
  lab_x <- ifelse(lab_x == "P_0", "P[0]", lab_x)
  
  all_df %>%
    filter(type == "result") %>%
    ggplot(aes(x = !!ensym(var_x), y = loglik))+
    geom_point(colour = point_colour)+
    theme_pubr() +
    geom_line(stat = "smooth", method = "loess", formula = y ~ x, alpha = 0.25,
              colour = point_colour, span = 0.75, size = 1) +
    labs(
      x = parse(text = lab_x))
}

# Axis text size (ats)
profile_plot <- function(profile_df, prof_var, maxloglik, ci_cutoff, est_colour,
                         ats = 11) {
  
  prof_var <- as.character(prof_var)
  
  profile_df %>%
    filter(is.finite(loglik)) %>%
    group_by(round(!!ensym(prof_var),5)) %>%
    filter(rank(-loglik)< 3) %>%
    ungroup() %>%
    ggplot(aes(x = !!ensym(prof_var),y = loglik))+
    geom_point(colour = est_colour, size = 0.5)+
    geom_line(stat = "smooth", method = "loess", formula = y ~ x, alpha = 0.25,
              colour = est_colour, span = 0.75, size = 0.8) +
    geom_hline(color = "red", yintercept = ci_cutoff, linetype = "dashed")+
    lims(y = maxloglik - c(5,0))+
    labs(y = "Log-lik", x = parse(text = prof_var),
         subtitle = "Unadjusted profile") +
    theme_pubr() +
    theme(axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60", size = ats),
          axis.ticks = element_line(colour = "grey60"))
}

plot_MCAP <- function(prof_ll, var_name, span = 0.75, est_colour,
                      axis_text_size = 11){
  
  var_name <- as.character(var_name)
  
  ggplot(prof_ll, aes(x = !!ensym(var_name), y = loglik)) +
    geom_point(alpha = 0.7, colour = est_colour, size = 0.3) +
    geom_line(stat = "smooth", method = "loess", formula = y ~ x, alpha = 0.25,
              colour = est_colour, span = span, size = 0.5) +
    geom_hline(yintercept = cut_off_smth, colour = "red", linetype = "dotted") +
    labs(y = "Log-lik", x = parse(text = var_name),
         subtitle = "MCAP") +
    theme_pubr() +
    theme(axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60", size = axis_text_size),
          axis.ticks = element_line(colour = "grey60"))
}

plot_guesses <- function(guesses_df, point_size = 1, point_colour) {
  aes_points <- list(continuous = wrap("points", alpha = 0.5, size = point_size,
                                       colour = point_colour))
  
  guesses_df <- rename(guesses_df, `P[0]` = `P_0`)
  
  ggpairs(guesses_df, 
          upper = aes_points,
          lower = aes_points,
          diag  = list(continuous =  'blankDiag'),
          labeller = label_parsed) +
    theme_pubr() +
    theme(axis.text = element_text(size = 6))
}

plot_priors <- function(){
  g1 <- ggplot(NULL, aes(c(0, 10))) + 
    geom_area(stat = "function", fun = dlnorm, fill = "grey95", 
              colour = "grey60") +
    scale_x_continuous(breaks = c(0, 5, 10)) +
    theme_pubr() +
    labs(y        = "Probability density",
         x        = bquote(zeta),
         subtitle = "lognormal(0,1)") +
    theme(axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 8),
          plot.subtitle = element_text(size = 8, colour = "grey60"))
  
  g2 <- g1 +
    labs(x = "P(0)", y = "", subtitle = "lognormal(0,1)")
  
  g3   <- ggplot(NULL, aes(c(0, 1))) + 
    geom_area(stat = "function", fun = dbeta, fill = "grey95", 
              colour = "grey60", args = list(shape1 = 2, shape2 = 2)) +  
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    theme_pubr() +
    labs(y        = "Probability density",
         x        = bquote(nu),
         subtitle = "beta(2,2)") +
    theme(axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 8),
          plot.subtitle = element_text(size = 8, colour = "grey60"))
  
  g4 <- g3 + labs(x = bquote(upsilon), y = "", subtitle = "beta(2,2)")
  
  (g1 + g2) / (g3 + g4)
}

plot_lik_surface <- function(tidy_ll_df, point_colour) {
  
  tidy_ll_df <- mutate(tidy_ll_df, 
                       name = ifelse(name == "P_0","P[0]", name))
  
  ggplot(tidy_ll_df, aes(x = value, y = loglik)) +
    geom_point(colour = point_colour, alpha = 0.75) +
    facet_wrap(~name, scales = "free", labeller = label_parsed) +
    geom_hline(yintercept = cutoff, linetype = "dashed", colour = "red") +
    theme_pubr() 
}

plot_demo <- function(demo_df, title, colour) {
  
  non_h   <- demo_df %>% filter(highlight == FALSE)
  high_df <- demo_df %>% filter(highlight == TRUE)
  ggplot(non_h, aes(x = week, y = beta)) +
    geom_line(aes(group = .id), colour = colour, alpha = 0.1) +
    geom_line(data = high_df, aes(group = .id), colour = colour) +
    scale_y_continuous(limits = c(0, 10)) +
    labs(title = title,
         y     = parse(text = "beta[t]"),
         x     = "Week") +
    theme_pubr() +
    theme(axis.line = element_line(colour = "grey65"),
          axis.ticks = element_line(colour = "grey65"),
          axis.text = element_text(colour = "grey65"))
}

plot_time_comparison <- function(time_df, tt){
  ggplot(time_df, aes(x = as.factor(order), y = time)) +
    geom_lollipop(point.colour = STH_colour, colour = "grey70") +
    coord_flip() +
    theme_pubr() +
    labs(y = "Elapsed time [mins]", x = "Delay order",
         caption = str_glue("Total computational time: {tt} mins")) +
    theme(axis.line = element_line(colour = "grey75"),
          axis.ticks = element_line(colour = "grey75"),
          axis.text = element_text(colour = "grey45"),
          axis.title = element_text(colour = "grey40"))
}


plot_traces <- function(sf) {
  traces_cols <- c("#D8DDB7", "#F5F3E7", "#1B4D60", "#578372")
  
  traceplot(sf, pars = c("zeta", "nu", "upsilon", "P_0")) +
    scale_colour_manual(values = traces_cols) +
    facet_wrap(~parameter, labeller = label_parsed, scales = "free")
}


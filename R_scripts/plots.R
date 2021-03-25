library(GGally)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(scales)

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
    geom_col(colour = "white", fill = "steelblue", alpha = 0.95) +
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
         title = "Daily epicurve",
         subtitle = str_glue("Total cases: {formatted_tc}")) +
    theme(plot.subtitle = element_text(size = 8, colour = "grey60"),
          plot.title = element_text(size = 9, colour = "grey39"),
          axis.title = element_text(size = 8, colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60", size = 6),
          axis.ticks = element_line(colour = "grey60"))
}

daily_epi_trend <- function(irish_data) {
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
         title = "Trend") +
    geom_point(colour = "steelblue") +
    stat_smooth(se = FALSE, span = 0.25, colour = "lightskyblue3") +
    theme_pubclean() +
    theme(axis.title = element_text(size = 8, colour = "grey40"),
          axis.text  = element_text(colour = "grey60", size = 6),
          plot.title = element_text(size = 9, colour = "grey39"),
          axis.ticks = element_line(colour = "grey60"))
}

weekly_epicurve <- function(wkl_data) {
  ggplot(wkl_data, aes(x = week, y = y)) +
    geom_col(fill = "steelblue") +
    scale_x_continuous(breaks = 0:12) +
    theme_pubclean() +
    labs(y = "Incidence [new cases per week]",
         x = "Week of Lab specimen collection",
         title = "Weekly epicurve") +
  theme_pubr() +
    theme(plot.title = element_text(size = 9, colour = "grey39"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          axis.title = element_text(size = 8, colour = "grey40"))
}

incidence_graphs <- function(irish_data, formatted_tc, wkl_inc) {
  g1 <- daily_epicurve(irish_data, formatted_tc)
  g2 <- daily_epi_trend(irish_data)
  g3 <- weekly_epicurve(wkl_inc)
  
  (g1/g2) | g3
}

driving_missing_data <- function(raw_data) {
  ggplot_na_distribution(raw_data) +
    labs(title = "Ireland's driving data",
         y = "Index",
         caption = "This dataset has been normalised by its initial value (141)") +
    theme_pubclean() +
    theme(plot.caption = element_text(colour = "grey70"))
}

imputed_driving_data <- function(raw_data, imputed_data) {
  ggplot_na_imputations(raw_data, imputed_data) +
    labs(title = "Ireland's driving data",
         y = "Index",
         caption = "Missing data has been replaced by linear interpolation estimates") +
    theme_pubclean() +
    theme(plot.caption = element_text(colour = "grey70"))
}

loglik_traces <- function(traces_df, lims) {
  traces_df %>% filter(variable == "loglik") %>% 
    ggplot(aes(x = iteration,y = value,group = L1))+
    geom_line(colour = "purple", alpha = 0.25)+
    scale_y_continuous(limits = lims) +
    guides(color=FALSE) +
    labs(x = "Iteration", y = "Log-likelihood") +
    theme_pubclean()
}

loglik_se_sens <- function(sens_results) {
  ggplot(sens_results, aes(x = Np, y = loglik.se)) +
    scale_x_log10(labels = comma) +
    geom_line(colour = "grey90") +
    geom_smooth(se = FALSE) +
    geom_point() +
    theme_pubclean() +
    labs(y = "Log-likelihood SE", x = "Number of particles (Log scale)")
}

mif_traces <- function(traces_df, unknown_pars) {
  traces_df %>% 
    filter(variable %in% c("loglik", unknown_pars)) %>%
    ggplot(aes(x= iteration,y = value, group = L1))+
    geom_line(color = "steelblue", alpha = 0.25) +
    theme_pubr() +
    scale_y_continuous(labels = comma) +
    labs(y = "Value", x = "Iteration") +
    facet_wrap(~variable,scales = "free_y", labeller = label_parsed)
}

raw_likelihood <- function(all_df, var_x) {
  all_df %>%
    filter(type == "result") %>%
    ggplot(aes(x = !!ensym(var_x), y = loglik))+
    geom_point()+
    theme_pubr() +
    geom_smooth(se = FALSE) +
    labs(
      x = parse(text = var_x),
      title = "Raw profile likelihood"
    )
}

profile_plot <- function(profile_df, prof_var, ci_cutoff) {
  profile_df %>%
    filter(is.finite(loglik)) %>%
    group_by(round(!!ensym(prof_var),5)) %>%
    filter(rank(-loglik)< 3) %>%
    ungroup() %>%
    ggplot(aes(x = !!ensym(prof_var),y = loglik))+
    geom_point()+
    geom_smooth(method = "loess", se = FALSE)+
    geom_hline(color = "red", yintercept = ci_cutoff, linetype = "dashed")+
    lims(y=maxloglik-c(5,0))+
    labs(y = "Log-likelihood", x = parse(text = prof_var)) +
    theme_pubr()
}



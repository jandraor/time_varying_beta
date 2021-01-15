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
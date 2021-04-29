get_unadjusted_profile <- function(df, var_name) {
  
  maxloglik <- max(df$loglik, na.rm = TRUE)
  cutoff    <- maxloglik - 0.5 * qchisq(df = 1,p = 0.95)
  
  mle <- filter(df, loglik == maxloglik) %>% pull(var_name)
  
  ci  <- filter(df, loglik >= cutoff) %>% pull(var_name) %>% range
  
  list(estimates = 
         list(ll = ci[[1]],
              ul = ci[[2]],
              mle = mle),
       cutoff    = cutoff,
       maxloglik = maxloglik)
}

collate_prof_estimates <- function(est_list) {
  imap_dfr(est_list, function(vals_list, name) {
    data.frame(name = name,
               mle  = vals_list$mle,
               ll = vals_list$ll,
               ul = vals_list$ul) 
  })
}

get_par_summary <- function(est_list) {
  source("./R_scripts/R_estimates.R")
  
  par_summary <- collate_prof_estimates(est_list)
  
  zeta_estimates <- par_summary %>% filter(name == "zeta") %>% 
    select(-name) %>% unlist()
  
  R_estimates <- estimate_r(zeta_estimates) %>% bind_rows()
  R_row       <- data.frame(name = "R(0)") %>% bind_cols(R_estimates)
  
  par_summary <- bind_rows(par_summary, R_row) %>% 
    mutate(name = ifelse(name == "P_0", "P(0)", name)) # formatting
  
  par_summary
}

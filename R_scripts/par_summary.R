
summarise_pars <- function(par_names, cutoff_list) {
  map_df(par_names, function(par_name) {
    cutoff_CI <- cutoff_list[[par_name]]
    
    profile_par <- prof_est_df %>% 
      filter(profile_over == par_name) %>% 
      filter(loglik > cutoff_CI)
    
    lims <- range(profile_par[[par_name]] %>% as.numeric)
    
    data.frame(Parameter = par_name,
               lower_limit = lims[[1]], upper_limit = lims[[2]])
  })
}

## Taken from 'add reference'
## Weighted quantile function
wquant <- function (x, weights, probs = c(0.025,0.25, 0.5, 0.75, 0.975)) {
  idx <- order(x)
  x <- x[idx]
  weights <- weights[idx]
  w <- cumsum(weights)/sum(weights)
  rval <- approx(w,x,probs,rule=1)
  rval$y
}
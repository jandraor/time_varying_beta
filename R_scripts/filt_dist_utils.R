

summarise_filt_distr <- function(filt_states, var_name) {
  
  imap_dfr(filt_states, function(filt_dist_t, i, var_name) {
    vals <-filt_dist_t[var_name, ]
    qs   <- quantile(vals, c(0.025, 0.25, 0.5, 0.75, 0.975))
    data.frame(time = as.numeric(i), q2.5 = qs[[1]], q25 = qs[[2]], q50 = qs[[3]], 
               q75 = qs[[4]], q97.5 = qs[[5]])
  }, var_name = var_name)
}



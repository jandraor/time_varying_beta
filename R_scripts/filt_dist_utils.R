# Sample from the filtering distribution
sample_filt_dist <- function(var_input, hindcast_df, seed_val) {
  
  set.seed(seed_val)
  
  var_df <- hindcast_df %>% filter(var == var_input)
  
  map_df(seq(7, 77, 7), function(i){
    time_df <- filter(var_df, time == i) %>% 
      mutate(rel_weight = weight / sum(weight))
    
    samples <- sample(time_df$value, 1e4, replace = TRUE,
                      prob = time_df$rel_weight)
    gc()
    
    data.frame(time = i, value = samples)
  }) %>% mutate(week = time / 7)
}

summarise_filt_distr <- function(filt_states, var_name) {
  
  imap_dfr(filt_states, function(filt_dist_t, i, var_name) {
    vals <-filt_dist_t[var_name, ]
    qs   <- quantile(vals, c(0.025, 0.25, 0.5, 0.75, 0.975))
    data.frame(time = as.numeric(i), q2.5 = qs[[1]], q25 = qs[[2]], q50 = qs[[3]], 
               q75 = qs[[4]], q97.5 = qs[[5]])
  }, var_name = var_name)
}



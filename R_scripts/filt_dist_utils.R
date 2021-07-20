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

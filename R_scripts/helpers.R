
create_2d_slices <- function(x1, x2, params) {
  
  params_list <- as.list(params)
  
  params_list[[x1$name]] <-  rep(seq(from = x1$min, to = x1$max, length = 40), 
                                 each = 3)
  
  params_list[[x2$name]] <-  rep(seq(from = x2$min, to = x2$max, length = 40), 
                                 each = 3)
  
  expand.grid(params_list)
}

discrete_net_change <- function(sim_df, cumulative_var) {
   temp_df       <- sim_df[, c("time", cumulative_var)]
   temp_df       <- filter(temp_df, time - trunc(time) == 0)
   cml_vals      <- temp_df[ , cumulative_var]
   temp_df$value <- round(cml_vals  - lag(cml_vals ), 0)
   temp_df       <- slice(temp_df, -1)
   
   temp_df[, c("time", "value")]
}

calculate_time <- function(t_list) {
  t_obj <- t_list[[1]]
  (t_obj$toc - t_obj$tic) / 60
}
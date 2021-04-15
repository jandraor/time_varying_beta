
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

extract_mif_results <- function(mif_output) {
  n_iter <- length(mif_output$mif_results)
  
  sapply(mif_output$mif_results, function(x) class(x))-> class_vector
  
  valid_indexes <- which(class_vector == "mif2d_pomp")
  
  valid_results <- mif_output$mif_results[valid_indexes]
  
  n_valid <- length(valid_results)
  
  message(str_glue("{n_valid} out of {n_iter} are valid"))
  
  # needs to be improved
  
  mif_output <- c(valid_results[[1]], valid_results[[2]])
  
  for(i in 3:length(valid_results)) {
    mif_output <- c(mif_output, valid_results[[i]])
  }
  
  mif_output
}

extract_ll_df <- function(ll_output) {
  n_iter <- length(ll_output$ll_results)
  
  map_lgl(ll_output$ll_results, function(x) {
    result_class <- class(x)
    
    "data.frame" %in% result_class
  }) -> valid_indexes
  
  valid_results <- ll_output$ll_results[valid_indexes]
  
  n_valid <- length(valid_results)
  
  message(str_glue("{n_valid} out of {n_iter} are valid"))
  
  do.call("rbind", valid_results)
}

create_stan_file <- function(stan_text, filename) {
  file_connection <- file(filename)
  writeLines(stan_text, file_connection)
  close(file_connection)
}
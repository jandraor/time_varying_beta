create_stan_file <- function(stan_text, filename) {
  file_connection <- file(filename)
  writeLines(stan_text, file_connection)
  close(file_connection)
}

stock_inits <- function(mdl, manual_values) {
  stocks <- sd_stocks(mdl) %>% 
    mutate(init_value = as.character(init_value))
  
  for(i in 1:length(manual_values)) {
    val_list <- manual_values[[i]]
    
    j            <- which(stocks$name == val_list$name)
    stocks[j, 2] <- val_list$value
  }
  
  stock_inits <- str_glue("  y0[{1:nrow(stocks)}] = {stocks$init_value};") %>% 
    paste(collapse = "\n")
}

run_stan_file <- function(fn, fit_options, stan_fn) {
  
  if(!file.exists(fn)) {
    
  
    mod <- cmdstan_model(stan_fn)
    
    fit <- mod$sample(data            = fit_options$stan_d,
                      seed            = fit_options$seed,
                      chains          = 4,
                      parallel_chains = 4,
                      iter_warmup     = fit_options$warmup,
                      iter_sampling   = fit_options$sampling,
                      refresh         = 5,
                      save_warmup     = TRUE,
                      step_size       = fit_options$step_size,
                      adapt_delta     = fit_options$adapt_delta)  
      
  sf  <- rstan::read_stan_csv(fit$output_files())
  
  saveRDS(sf, fn)
    
  } else {
    sf <- readRDS(fn)
  }
  sf
}

construct_incidence_df <- function(posterior_df) {
  extract_timeseries_var("y1_hat", posterior_df) %>% 
    mutate(order = dly_o)
}

mase_per_iter <- function(sim_df, i, data_vector) {
  sim_df %>% group_by(iter) %>% 
    summarise(mase = mase(data_vector, value)) %>% 
    mutate(order = i)
}
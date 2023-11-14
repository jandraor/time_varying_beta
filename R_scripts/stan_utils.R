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


#' Run Stan file
#'
#' @param fn A string that indicates the path where file will be saved
#' @param fit_options A list
#' @param stan_fn A string that indicates Stan file's path
run_stan_file <- function(order, fit_options, backup_folder, stan_folder) {
  
  fn <- file.path(backup_folder, str_glue("stan_fit_order_{order}.rds"))
  
  if(!file.exists(fn)) {
    
    stan_fn <- file.path(stan_folder, str_glue("SEI3R_{order}_smth.stan"))
    
    tic.clearlog()
    tic()
  
    mod <- cmdstan_model(stan_fn)
    
    n_chains <- 4
    
    if(!is.null(fit_options$chains)) n_chains <- fit_options$chains
    
    print(fit_options$init)
    
    fit <- mod$sample(data            = fit_options$stan_d,
                      seed            = fit_options$seed,
                      chains          = n_chains,
                      parallel_chains = 4,
                      iter_warmup     = fit_options$warmup,
                      iter_sampling   = fit_options$sampling,
                      refresh         = 5,
                      save_warmup     = FALSE,
                      step_size       = fit_options$step_size,
                      adapt_delta     = fit_options$adapt_delta,
                      init            = fit_options$init)  
    
    toc(quiet = FALSE, log = TRUE)
    log.lst <- tic.log(format = FALSE)
    
    diag_path <- file.path(backup_folder, str_glue("diag_{order}.txt"))
    diagnosis <- fit$cmdstan_diagnose()
    writeLines(diagnosis$stdout, diag_path)
    
    sf  <- rstan::read_stan_csv(fit$output_files())
    
    results <- list(sf = sf, time = log.lst) 
  
    saveRDS(results, fn)
    
  } else {
    results <- readRDS(fn)
  }
  results
}

run_stan_files <- function(n_orders, stan_d, folder, ad_list, sz_list, mt_list) {
  
  posterior_list <- lapply(1:n_orders, \(dly_o) {
    
    ad <- ad_list[[dly_o]]
    sz <- sz_list[[dly_o]]
    mt <- mt_list[[dly_o]]
    
    fn <- file.path(folder, str_glue("stan_fit_order_{dly_o}.rds"))
    
    if(!file.exists(fn)) {
      
      stan_path <- str_glue("./Stan_files/SEI3R_{dly_o}_smth.stan")
      mod       <- cmdstan_model(stan_path)
      
      message(str_glue("Delay: {dly_o}"))
      message(str_glue("Step size: {sz}"))
      message(str_glue("Adapt delta: {ad}"))
      message(str_glue("Max treedepth: {mt}"))
      
      tic.clearlog()
      tic()
      
      fit <- mod$sample(data            = stan_d,
                        chains          = 4,
                        parallel_chains = 4,
                        iter_warmup     = 2000,
                        iter_sampling   = 2000,
                        refresh         = 100,
                        save_warmup     = FALSE,
                        adapt_delta     = ad,
                        step_size       = sz,
                        max_treedepth   = mt)
      
      toc(quiet = FALSE, log = TRUE)
      log.lst <- tic.log(format = FALSE)
      
      diag_path <- file.path(backup_folder, str_glue("diag_{dly_o}.txt"))
      diagnosis <- fit$cmdstan_diagnose()
      writeLines(diagnosis$stdout, diag_path)
      
      posterior_df <- as_draws_df(fit$draws())
      
      results <- list(sf = posterior_df, time = log.lst) 
      
      saveRDS(results, fn)
      
    } else {
      results <- readRDS(fn)
    }
    
    results
  })
}

construct_incidence_df <- function(posterior_df, dly_o) {
  extract_timeseries_var("delta_x_1", posterior_df) %>% 
    mutate(order = dly_o)
}

mase_per_iter <- function(sim_df, i, data_vector) {
  sim_df %>% group_by(iter) %>% 
    summarise(mase = mase(data_vector, value)) %>% 
    mutate(order = i)
}
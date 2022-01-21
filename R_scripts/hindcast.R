hindcast <- function(unk_pars_df, pomp_mdl, fixed_pars, folder, seed, 
                     n_cores) {
  
  hc_vars <- c("C", "Z", "Re")
  
  sapply(hc_vars, function(hc_var) {
    filename <- str_glue("{folder}/{hc_var}.rds")
    file.exists(filename)
  }) |> sum() -> flnm_vals 
  
  time_file <- str_glue("{folder}/hindcast_time.rds")
  
  
  if(flnm_vals != 3) {
    
    registerDoRNG(seed)
    registerDoParallel(cores = n_cores)
    
    message(paste0("Number of working cores: "), getDoParWorkers())
    
    message(paste0("Starting hindcast at ", date()))
    
    tic.clearlog()
    tic()
    
    foreach(unk_pars = iter(unk_pars_df,"row"), .combine = c) %dopar% {
      
      library(dplyr)
      library(pomp)
      library(purrr)
      library(reshape2)
      library(tibble)
      library(tidyr)
      
      source("./R_scripts/R_estimates.R")
      
      pf <- pfilter(pomp_mdl, params = c(fixed_pars, unlist(unk_pars)),
                    Np = 1e4, save.states = TRUE)
      
      ss <- saved.states(pf)
      
      id <- as.numeric(rownames(unk_pars))
      
      imap_dfr(ss, function(states, i) {
        
        trans_states_df <- states[c("C", "Z", "S"), ] |>  t() |> 
          as.data.frame()
        
        calculations_df <- trans_states_df |> 
          mutate(s = S / fixed_pars[["N"]],
                 zeta = unk_pars[["zeta"]],
                 R  = estimate_r(zeta * Z),
                 Re = R * s) |>  
          select(C, Z, R, Re) |>  
          mutate(rep = row_number())
        
        tidy_df <- calculations_df |> 
          pivot_longer(-rep, names_to = "var", values_to = "value")
        
        tidy_df <- tidy_df |>  
          mutate(time   = as.numeric(i),
                 id     = id,
                 loglik = logLik(pf))
      }) -> iter_results
      
      saveRDS(iter_results, stringr::str_glue("{folder}/hc_iter_{id}.rds"))
      
      1
    } -> indicator
    
    n_iters <- nrow(unk_pars_df)
    
    walk(hc_vars, function(hc_var) {
      
      print(hc_var)
      
      var_df <- map_df(1:n_iters, function(id) {
        
        iter_file <- str_glue("{folder}/hc_iter_{id}.rds")
        df        <- readRDS(iter_file) |> filter(var == hc_var)
      })
      
      saveRDS(var_df, str_glue("{folder}/{hc_var}.rds"))
    })
    
    walk(1:n_iters, function(id) unlink(str_glue("{folder}/hc_iter_{id}.rds")))
    
    toc(quiet = FALSE, log = TRUE)
    log.lst   <- tic.log(format = FALSE)
    time_list <- list(time = log.lst)
    saveRDS(time_list, time_file)
  } else {
    time_list <- readRDS(time_file)
  }
  
  time_list
}
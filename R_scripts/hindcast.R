hindcast <- function(unk_pars_df, pomp_mdl, fixed_pars, filename, seed, 
                     n_cores) {
  
  if(!file.exists(filename)) {
    registerDoRNG(seed)
    registerDoParallel(cores = n_cores)
    
    message(paste0("Number of working cores: "), getDoParWorkers())
    
    message(paste0("Starting hindcast at ", date()))
    
    tic.clearlog()
    tic()
    
    foreach(unk_pars = iter(unk_pars_df,"row"), .combine = rbind) %dopar% {
      
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
      
      imap_dfr(ss, function(states, i) {
        trans_states_df <- states[c("C", "Z", "S"), ] %>% t() %>% 
          as.data.frame()
        
        calculations_df <- trans_states_df %>% 
          mutate(s = S / fixed_pars[["N"]],
                 zeta = unk_pars[["zeta"]],
                 R  = estimate_r(zeta * Z),
                 Re = R * s) %>% 
          select(C, Z, R, Re) %>% 
          mutate(rep = row_number())
        
        tidy_df <- calculations_df %>% 
          pivot_longer(-rep, names_to = "var", values_to = "value")
        
        tidy_df <- tidy_df %>% 
          mutate(time = as.numeric(i),
                 id     = as.numeric(rownames(unk_pars)),
                 loglik = logLik(pf))
      })
    } -> hc_results
    
    toc(quiet = FALSE, log = TRUE)
    log.lst <- tic.log(format = FALSE)
    results <- list(hindcast = hc_results, time = log.lst)
    saveRDS(results, filename)
  } else {
    results <- readRDS(filename)
  }
  
  results
}
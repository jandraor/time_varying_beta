global_search <- function(guesses, fixed_params, mf1, fn, seed, n_cores) {
  
  if(!file.exists(fn)) {
    
    registerDoParallel(cores = n_cores)
    registerDoRNG(seed)
    
    tic.clearlog()
    tic()
    
    foreach(guess=iter(guesses,"row"), .combine = c) %dopar% {
      mf1 %>%
        mif2(params = c(unlist(guess), fixed_params)) %>%
        mif2(Nmif = 100) -> mf
    } -> mf_results
    
    toc(quiet = FALSE, log = TRUE)
    log.lst <- tic.log(format = FALSE)
    results <- list(mf_results = mf_results, time = log.lst)
    saveRDS(results, fn)
  } else {
    results <- readRDS(fn)  
  }
  
  results
}


# replicate(
#   10,
#   mf %>% pfilter(Np = 100000) %>% logLik()
# ) %>%
#   logmeanexp(se=TRUE) -> ll
# 
# mf %>% coef() %>% bind_rows() %>%
#   bind_cols(loglik = ll[1],loglik.se=ll[2])
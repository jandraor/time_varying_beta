global_search <- function(guesses, fixed_params, mf1, fn, seed, n_cores) {
  
  if(!file.exists(fn)) {
    message("...Starting global search...")
    
    registerDoParallel(cores = n_cores)
    registerDoRNG(seed)
    
    message(paste0("Number of working cores: "), getDoParWorkers())
    
    tic.clearlog()
    tic()
    
    foreach(guess = iter(guesses,"row"), .combine = c, 
            .errorhandling = 'pass') %dopar% {
      library(dplyr)
      library(pomp)
      
      out <- tryCatch({
        mf1 %>%
          mif2(params = c(unlist(guess), fixed_params)) %>%
          mif2(Nmif = 100) -> mf
        
        list(mf)
      },
      error = function(cond) {
        error_list <- list(list(guess = unlist(guess),
                           error = cond))
        return(error_list)
      })
      
      return(out)
    } -> mf_results
    
    toc(quiet = FALSE, log = TRUE)
    log.lst <- tic.log(format = FALSE)
    results <- list(mif_results = mf_results, time = log.lst)
    saveRDS(results, fn)
  } else {
    results <- readRDS(fn)  
  }
  
  results
}
estimate_avg_lik <- function(pomp_mdl, seed, cores, Np = 5000) {
  registerDoParallel(cores = 6)
  registerDoRNG(652643293)
  
  foreach (i = 1:10, .combine = c) %dopar% {
    library(dplyr)
    library(pomp)
    pomp_mdl %>% pfilter(Np = Np)
  } -> pf
  logLik(pf) -> ll
  logmeanexp(ll, se=TRUE)
}

iter_filt_profile <- function(mf1, guesses, fixed_params, perturbations, 
                              filename, seed, n_cores = 7) {
  
  if(!file.exists(filename)) {
    
    registerDoRNG(seed)
    registerDoParallel(cores = n_cores)
    
    tic.clearlog()
    tic()
    
    foreach(guess = iter(guesses,"row"), .combine = c) %dopar% {
      library(dplyr)
      library(pomp)
      
      mf1 %>%
        mif2(params = c(unlist(guess),fixed_params),
             rw.sd = perturbations) %>%
        mif2(Nmif = 100, cooling.fraction.50 = 0.3) -> mf
      

    } -> mif_results
    toc(quiet = FALSE, log = TRUE)
    log.lst <- tic.log(format = FALSE)
    results <- list(mif_results = mif_results, time = log.lst)
    saveRDS(results, filename)
  } else {
    results <- readRDS(filename)
  }
  results
}

mif_ll <- function(mf_list, Np = 100000, seed, n_cores, filename) {
  
  if(!file.exists(filename)) {
    
    registerDoRNG(seed)
    registerDoParallel(cores = n_cores)
    
    message(paste0("Number of working cores: "), getDoParWorkers())
    
    tic.clearlog()
    tic()

    foreach(mf = mf_list, .combine = rbind) %dopar% {
      
      library(dplyr)
      library(pomp)
      
      replicate(10, mf %>% pfilter(Np = Np) %>% logLik()) %>%
      logmeanexp(se = TRUE) -> ll
      
      mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik = ll[1],loglik.se = ll[2])
    } -> ll_results
    
    toc(quiet = FALSE, log = TRUE)
    log.lst <- tic.log(format = FALSE)
    results <- list(ll_results = ll_results, time = log.lst)
    saveRDS(results, filename)
  } else {
    results <- readRDS(filename)
  }
}



estimate_avg_lik <- function(pomp_mdl, seed, cores, Np = 5000,
                             n_iter = 10) {
  registerDoParallel(cores = cores)
  registerDoRNG(seed)
  
  foreach (i = 1:n_iter, .combine = c) %dopar% {
    library(dplyr)
    library(pomp)
    pomp_mdl %>% pfilter(Np = Np)
  } -> pf
  logLik(pf) -> ll
  logmeanexp(ll, se = TRUE)
}

iter_filt_profile <- function(mf1, guesses, fixed_params, perturbations, 
                              filename, seed, n_cores = 7) {
  
  if(!file.exists(filename)) {
    
    registerDoRNG(seed)
    registerDoParallel(cores = n_cores)
    
    message(paste0("Number of working cores: "), getDoParWorkers())
    
    tic.clearlog()
    tic()
    
    foreach(guess = iter(guesses,"row"), .combine = c,
            .errorhandling = 'pass') %dopar% {
      library(dplyr)
      library(pomp)
              
      out <- tryCatch({
        mf1 %>%
          mif2(params = c(unlist(guess),fixed_params),
               rw.sd = perturbations) %>%
          mif2(Nmif = 100, cooling.fraction.50 = 0.3) -> mf
        list(mf)
      },
      error = function(cond) {
        error_list <- list(list(guess = unlist(guess),
                                error = cond))
        return(error_list)
      })
      
      return(out)
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
    
    message(paste0("Starting at ", date()))
    
    tic.clearlog()
    tic()

    foreach(mf = mf_list, .combine = c,
            .errorhandling = 'pass') %dopar% {
      
      library(dplyr)
      library(pomp)
              
      out <- tryCatch({
        replicate(10, mf %>% pfilter(Np = Np) %>% logLik()) %>%
          logmeanexp(se = TRUE) -> ll
        
        mf %>% coef() %>% bind_rows() %>%
          bind_cols(loglik = ll[1],loglik.se = ll[2]) -> ll_df
        list(ll_df)
      },
      error = function(cond) {
        error_list <- list(list(error = cond))
        return(error_list)
      })
    
    } -> ll_results
    
    toc(quiet = FALSE, log = TRUE)
    log.lst <- tic.log(format = FALSE)
    results <- list(ll_results = ll_results, time = log.lst)
    saveRDS(results, filename)
  } else {
    results <- readRDS(filename)
  }
  results
}

pf_sensitivity <- function(n_particles, pomp_mdl, n_cores, seed, fn,
                           n_iter) {
    
  if(!file.exists(fn)) {
  
  map_df(n_particles, function(Np) {
    tic.clearlog()
    tic()
    lik_result <- estimate_avg_lik(pomp_mdl, seed, n_cores, Np, n_iter)
    toc(quiet = FALSE, log = TRUE)
    log.lst <- tic.log(format = FALSE)
    
    data.frame(Np = Np, loglik = lik_result[[1]], loglik.se = lik_result[[2]],
               time = calculate_time(log.lst))
  }) -> results
  
  saveRDS(results, fn)
    
  } else {
    results <- readRDS(fn)
  }
  results
}

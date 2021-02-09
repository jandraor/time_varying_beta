likelihood_slice <- function(combinations, folder) {
  ll_slice_list <- purrr::map(combinations, function(indexes) {
    i       <- indexes[[1]]
    j       <- indexes[[2]]
    x1_name <- params_bounds[[i]]$name
    x2_name <- params_bounds[[j]]$name
    
    file_path <- str_glue("{folder}surface_exploration_{x1_name}_{x2_name}.rds")
    
    if(!file.exists(file_path)) {
      p <- create_2d_slices(params_bounds[[i]], params_bounds[[j]],
                            params)    
      
      registerDoParallel(cores = detectCores() - 2)
      registerDoRNG(421776444)
      
      tic.clearlog()
      tic()
      
      foreach (theta = iter(p,"row"), .combine = rbind, .inorder=FALSE) %dopar%
        {
          SEI3R_GBM %>% pfilter(params = theta, Np = 5000) -> pf
          theta$loglik <- logLik(pf)
          theta
        } -> p
      
      toc(quiet = TRUE, log = TRUE)
      log.lst <- tic.log(format = FALSE)
      
      p_list  <- list(p = p, time = log.lst)
      saveRDS(p_list, file_path)
    } else {
      p_list <- readRDS(file_path)
    }
    p_list
  })
}


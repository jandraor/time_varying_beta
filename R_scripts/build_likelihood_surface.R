build_likelihood_surface <- function(folder, unk_pars) {
  
  ls_file <- "local_search_ll.rds"
  gs_file <- "Global_search_ll.rds"
  
  profile_files <- str_glue("ifp_{unk_pars}_ll.rds")
  file_paths    <- file.path(folder, c(ls_file, gs_file, profile_files))
  
  map_df(file_paths, function(file) {
    obj <- readRDS(file)
    
    ll_results <- extract_ll_df(obj) |> 
      select(c(unk_pars, "loglik", "loglik.se"))
  }) |> 
    filter(loglik.se < 1) |> select(-loglik.se)
}
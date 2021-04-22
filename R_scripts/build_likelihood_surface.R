build_likelihood_surface <- function(mdl, unk_pars) {
  
  folder  <- str_glue("./Saved_objects/Irish_data/SEI3R_{mdl}/weekly/mdl_2")
  
  ls_file <- "local_search_mdl2_ll.rds"
  gs_file       <- "Global_search_mdl2_ll.rds"
  
  profile_files <- str_glue("ifp_{unk_pars}_ll_mdl2.rds")
  file_paths    <- file.path(folder, c(ls_file, gs_file, profile_files))
  
  map_df(file_paths, function(file) {
    obj <- readRDS(file)
    
    ll_results <- extract_ll_df(obj) %>% 
      select(c(unk_pars, "loglik", "loglik.se"))
  }) %>% 
    filter(loglik.se < 1) %>% select(-loglik.se)
}
source("./R_scripts/model_setup_GBM.R")

var_name <- "tau"

fn_ifp  <- file.path(folder, str_glue("ifp_{var_name}.rds"))
ifp_obj <- readRDS(fn_ifp)

source("./R_scripts/helpers.R")
profile_mif_results <- extract_mif_results(ifp_obj)

source("./R_scripts/likelihood_funs.R")
fn              <- file.path(folder, str_glue("ifp_{var_name}_ll_mdl2.rds"))
profile_results <- mif_ll(profile_mif_results, seed = 520044828, 
                          n_cores = 7, filename = fn)
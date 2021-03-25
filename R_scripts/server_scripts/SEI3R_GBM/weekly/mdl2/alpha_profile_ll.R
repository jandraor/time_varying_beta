source("./R_scripts/model_setup_GBM.R")

fn_ifp  <- file.path(folder, "ifp_alpha.rds")
ifp_obj <- readRDS(fn_ifp)

source("./R_scripts/helpers.R")
alpha_profile_mif_results <- extract_mif_results(ifp_obj)

source("./R_scripts/likelihood_funs.R")
fn              <- file.path(folder, "ifp_alpha_ll_mdl2.rds")
profile_results <- mif_ll(alpha_profile_mif_results, seed = 872114710, 
                          n_cores = 7, filename = fn)
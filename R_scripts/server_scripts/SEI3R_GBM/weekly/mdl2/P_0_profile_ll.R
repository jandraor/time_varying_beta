source("./R_scripts/model_setup_GBM.R")

fn_ifp  <- file.path(folder, "ifp_P_0.rds")
ifp_obj <- readRDS(fn_ifp)

source("./R_scripts/helpers.R")
P_0_profile_mif_results <- extract_mif_results(ifp_obj)

source("./R_scripts/likelihood_funs.R")
fn              <- file.path(folder, "ifp_P_0_ll_mdl2.rds")
profile_results <- mif_ll(P_0_profile_mif_results, seed = 990258595, 
                          n_cores = 7, filename = fn)
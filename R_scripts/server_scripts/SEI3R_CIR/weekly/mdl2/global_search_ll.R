source("./R_scripts/model_setup_CIR.R")

gs_obj      <- readRDS(file.path(folder, "Global_search_mdl2.rds"))
mifs_global <- extract_mif_results(gs_obj)

source("./R_scripts/likelihood_funs.R")
fn          <- file.path(folder, "Global_search_mdl2_ll.rds" )
ll_obj      <- mif_ll(mifs_global, Np = 100000, 1270401374, 
                      n_cores = 7, fn)
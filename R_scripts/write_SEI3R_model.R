write_SEI3R_model <- function(mdl_path, stan_path) {
  source("./R_scripts/stan_utils.R")
  
  mdl      <- read_xmile(mdl_path)
  stocks   <- sd_stocks(mdl)
  
  ODE_fn    <- "SEI3R"
  stan_fun  <- stan_ode_function(mdl_path, ODE_fn, 
                                 pars = c("zeta", "inv_nu", "upsilon"))
  
  fun_exe_line <- str_glue("  o = ode_rk45({ODE_fn}, y0, t0, ts, params);") 
  
  stan_data <- stan_data(c("y1"), type = c("int"), inits = FALSE)
  
  stan_params <- paste(
    "parameters {",
    "  real<lower = 0>            zeta;",
    "  real<lower = 0, upper = 1> nu;",
    "  real<lower = 0, upper = 1> upsilon;",
    "  real<lower = 0>            P_0;",
    "}", sep = "\n")
  
  manual_values <- list(list(name = "S", value = "4937794 - P_0"),
                        list(name = "P", value = "P_0"))
  
  inits <- stock_inits(mdl, manual_values)
  
  k <- which(stocks$name == "C")
  
  stan_tp <- paste(
    "transformed parameters{",
    "  vector[n_difeq] o[n_obs]; // Output from the ODE solver",
    "  real y1_hat[n_obs];",
    "  vector[n_difeq] y0;",
    "  real params[n_params];",
    inits,
    "  params[1] = zeta;",
    "  params[2] = 1 / nu;",
    "  params[3] = upsilon;",
    fun_exe_line,
    str_glue("  y1_hat[1] =  o[1, {k}]  - y0[{k}];"),
    "  for (i in 1:n_obs-1) {",
    str_glue("    y1_hat[i + 1] = o[i + 1, {k}] - o[i, {k}] + 1e-4;"),
    "  }",
    "}", sep = "\n")
  
  stan_model <- paste(
    "model {",
    "  zeta    ~ lognormal(0, 1);",
    "  upsilon ~ beta(2, 2);",
    "  nu      ~ beta(2, 2);",
    "  P_0     ~ lognormal(0, 1);",
    "  y1      ~ poisson(y1_hat);",
    "}",
    sep = "\n")
  
  stan_gc <- paste(
    "generated quantities {",
    "  real log_lik;",
    "  log_lik = poisson_lpmf(y1 | y1_hat);",
    "}", sep = "\n")
  
  stan_text   <- paste(stan_fun, stan_data, stan_params,
                       stan_tp, stan_model, stan_gc, sep = "\n")
  
  create_stan_file(stan_text, stan_path)
}

write_SEI3R_model2 <- function(mdl_path, stan_path) {
  source("./R_scripts/stan_utils.R")
  
  mdl      <- read_xmile(mdl_path)
  stocks   <- sd_stocks(mdl)
  
  ODE_fn    <- "SEI3R"
  stan_fun  <- stan_ode_function(mdl_path, ODE_fn, 
                                 pars = c("zeta", "inv_nu", "upsilon"))
  
  fun_exe_line <- str_glue("  o = ode_rk45({ODE_fn}, y0, t0, ts, params);") 
  
  stan_data <- stan_data(c("y1"), type = c("int"), inits = FALSE)
  
  stan_params <- paste(
    "parameters {",
    "  real<lower = 0>            zeta;",
    "  real<lower = 0, upper = 1> nu;",
    "  real<lower = 0, upper = 1> upsilon;",
    "  real<lower = 0>            P_0;",
    "  real<lower = 0>            phi;",
    "}", sep = "\n")
  
  manual_values <- list(list(name = "S", value = "4937794 - P_0"),
                        list(name = "P", value = "P_0"))
  
  inits <- stock_inits(mdl, manual_values)
  
  k <- which(stocks$name == "C")
  l <- which(stocks$name == "Z")
  
  stan_tp <- paste(
    "transformed parameters{",
    "  vector[n_difeq] o[n_obs]; // Output from the ODE solver",
    "  real y1_hat[n_obs];",
    "  vector[n_difeq] y0;",
    "  real params[n_params];",
    "  real inv_phi;",
    "  inv_phi = 1 / sqrt(phi);",
    inits,
    "  params[1] = zeta;",
    "  params[2] = 1 / nu;",
    "  params[3] = upsilon;",
    fun_exe_line,
    str_glue("  y1_hat[1] =  o[1, {k}]  - y0[{k}];"),
    "  for (i in 1:n_obs-1) {",
    str_glue("    y1_hat[i + 1] = o[i + 1, {k}] - o[i, {k}] + 1e-4;"),
    "  }",
    "}", sep = "\n")
  
  stan_model <- paste(
    "model {",
    "  zeta    ~ lognormal(0, 1);",
    "  upsilon ~ beta(2, 2);",
    "  nu      ~ beta(2, 2);",
    "  P_0     ~ lognormal(0, 1);",
    "  phi     ~ normal(0, 1);",
    "  y1      ~ neg_binomial_2(y1_hat, inv_phi);",
    "}",
    sep = "\n")
  
  stan_gc <- paste(
    "generated quantities {",
    "  real log_lik;",
    "  log_lik = neg_binomial_2_lpmf(y1 | y1_hat, inv_phi);",
    "}", sep = "\n")
  
  stan_text   <- paste(stan_fun, stan_data, stan_params,
                       stan_tp, stan_model, stan_gc, sep = "\n")
  
  create_stan_file(stan_text, stan_path)
}
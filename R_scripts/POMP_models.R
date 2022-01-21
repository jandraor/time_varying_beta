get_fixed_params <- function() {
  sigma_val   <- 1 / 3   # Incubation period
  eta_val     <- 1 / 2.1 # Duration of preclinical infectiousness
  gamma_val   <- 1 / 2.9 # Duration of clinical infectiousness
  kappa_val   <- 1 / 5 # Duration of subclinical infectiousness
  N_val       <- 4937796
  mu_val      <- 1 / 2 # Relative infectiousness
  omega_val   <- 0.73 # Clinical fraction
  
  c(N = N_val, mu = mu_val, omega = omega_val, gamma = gamma_val,
    sigma = sigma_val, eta = eta_val, kappa = kappa_val)
}

pomp_GBM <- function(incidence_mdl, mobility_mdl, obs_df, mdl_filename, dt) {
  
  mm_obj <- get_meas_mdl(incidence_mdl, mobility_mdl)
  
  rmeas <- mm_obj$rmeas
  dmeas <- mm_obj$dmeas
  
  Csnippet("
    S = N - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    Z = 1;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW     = rnorm(0,sqrt(dt));
    double lambda = zeta * Z * (I + P + mu * A) / N;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    Z+= alpha*Z*dW;
    C+= (eta*P)*dt;
  ") -> SEI3R_GBM_step
  
  
  par_obj   <- get_params("GBM", incidence_mdl, mobility_mdl)
  params    <- par_obj$all
  par_names <- names(params)
  
  o_names <- "y1"
  
  if(mobility_mdl) o_names <- c("y1", "y2")
  
  obs_df |>   
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_GBM_step, delta.t = dt),
      statenames = c("S", "E", "P","I", "R", "A", "C", "Z"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log = names(par_obj$unknown)),
      obsnames = o_names,
      cdir = ".", 
      cfile = mdl_filename) -> mdl
  
  list(mdl  = mdl,
       pars = par_obj)
}

pomp_CIR <- function(incidence_mdl, mobility_mdl, obs_df, mdl_filename, dt) {
  
  mm_obj <- get_meas_mdl(incidence_mdl, mobility_mdl)
  
  rmeas <- mm_obj$rmeas
  dmeas <- mm_obj$dmeas
  
  Csnippet("
    S = N - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    Z = 1;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW     = rnorm(0,sqrt(dt));
    double lambda = zeta * Z * (I + P + mu * A) / N;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    Z+= (nu * (upsilon - Z)) * dt + sqrt(alpha)*Z*dW;
    C+= (eta*P)*dt;
  ") -> SEI3R_CIR_step
  
  par_obj   <- get_params("CIR", incidence_mdl, mobility_mdl)
  params    <- par_obj$all
  par_names <- names(params)
  
  par_unk   <- names(par_obj$unknown)
  
  log_pars  <- par_unk[!par_unk %in% c("upsilon", "nu")]
  
  o_names <- "y1"
  
  if(mobility_mdl) o_names <- c("y1", "y2")
  
  obs_df |> 
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_CIR_step, delta.t = dt),
      statenames = c("S", "E", "P","I", "R", "A", "C", "Z"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log   = log_pars,
                                 logit = c("upsilon", "nu")),
      obsnames = o_names,
      cdir     = ".", 
      cfile    = mdl_filename
    ) -> mdl
  
  
  list(mdl  = mdl,
       pars = par_obj)
}

get_params <- function(PM, incidence_mdl, mobility_mdl) {
  
  zeta_guess    <- 3 / 5
  tau_guess     <- 1
  B_guess       <- zeta_guess
  P_guess       <- 1
  alpha_guess   <- 0.1
  nu_guess      <- 0.04
  upsilon_guess <- 0.1
  phi_guess     <- 0.3
  
  unknown_pars  <- c(zeta = zeta_guess, P_0 = P_guess, alpha = alpha_guess)
  
  if(PM == "CIR") {
    unknown_pars <- c(unknown_pars, c(nu = nu_guess, upsilon = upsilon_guess))
  }
  
  if(incidence_mdl == "Nbin") {
    unknown_pars <- c(unknown_pars, c(phi = phi_guess))
  }
  
  if(mobility_mdl) {
    unknown_pars <- c(unknown_pars, c(tau = tau_guess))
  }
  
  fixed_pars   <- get_fixed_params()
  
  params       <- c(fixed_pars, unknown_pars)
  
  list(fixed   = fixed_pars,
       unknown = unknown_pars,
       all     = params)
} 

get_meas_mdl <- function(incidence_mdl, mobility_mdl) {
  
  if(incidence_mdl == "Pois") rmeas_inc <- "y1  = rpois(C);"
  
  if(incidence_mdl == "Nbin") rmeas_inc <- "y1  = rnbinom_mu(C, 1 / phi);"
  
  if(!mobility_mdl) rmeas_text <- rmeas_inc
  
  if(mobility_mdl) {
    rmeas_text <- paste(rmeas_inc, "y2  = rnorm(Z, tau);", sep = "\n")
  }
  
  rmeas <- Csnippet(rmeas_text)
  
  #-----------------------------------------------------------------------------
  
  if(incidence_mdl == "Pois") dmeas_inc <- "lik = dpois(y1,C,give_log)"
  
  if(incidence_mdl == "Nbin") {
    dmeas_inc <- "lik = dnbinom_mu(y1, 1 / phi, C, give_log)"
  }
  
  if(!mobility_mdl) dmeas_text <- paste0(dmeas_inc, ";")
  
  if(mobility_mdl) {
    dmeas_text <- paste0(dmeas_inc, " + ", "dnorm(y2, Z, tau, give_log);")
  }
  
  dmeas <- Csnippet(dmeas_text)
  
  list(rmeas = rmeas,
       dmeas = dmeas)
}
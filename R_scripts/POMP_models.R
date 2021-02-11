get_POMP_model <- function(obs_df, params) {
  
  Csnippet("
    y = rpois(C);  
  ") -> rmeas
  
  Csnippet("
    lik = dpois(y,C,give_log) + 1e-5;
  ") -> dmeas
  
  Csnippet("
    S = 4999970 - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    B = B_0;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW     = rnorm(0,sqrt(dt));
    double lambda = B * (I + P + mu * A) / N;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    B+= alpha*B*dW;
    C+= (rho*eta*P)*dt;
  ") -> SEI3R_GBM_step
  
  par_names <- names(params)
  
  obs_df %>%  
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_GBM_step, delta.t = 0.01),
      statenames = c("S", "E", "P","I", "R", "A", "C", "B"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log   = c("B_0", "P_0"))
    ) -> SEI3R_GBM
  
}


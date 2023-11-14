functions {
  vector X_model(real time, vector y, array[] real params) {
    vector[8] dydt;
    real E_to_P;
    real E_to_A;
    real P_to_I;
    real A_to_R;
    real I_to_R;
    real lambda;
    real inv_nu;
    real P_to_C;
    real S_to_E;
    real adjust_Z;
    E_to_P = 0.73*0.3333333333*y[3];
    E_to_A = (1-0.73)*0.3333333333*y[3];
    P_to_I = 0.4761904762*y[4];
    A_to_R = 0.2*y[7];
    I_to_R = 0.3448275862*y[5];
    lambda = (params[1]*y[1]/4937796)*(y[5]+y[4]+0.5*y[7]);
    inv_nu = 1/0.08;
    P_to_C = P_to_I;
    S_to_E = lambda*y[2];
    adjust_Z = (0.05-y[1])/inv_nu;
    dydt[1] = adjust_Z;
    dydt[2] = -S_to_E;
    dydt[3] = S_to_E-E_to_P-E_to_A;
    dydt[4] = E_to_P-P_to_I;
    dydt[5] = P_to_I-I_to_R;
    dydt[6] = A_to_R+I_to_R;
    dydt[7] = E_to_A-A_to_R;
    dydt[8] = P_to_C;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  array[n_obs] int y1;
  real t0;
  array[n_obs] real ts;
  vector[8] x0;
}
parameters {
  real<lower = 0> zeta;
}
transformed parameters{
  array[n_obs] vector[8] x; // Output from the ODE solver
  array[1] real params;
  array[n_obs] real delta_x_1;
  params[1] = zeta;
  x = ode_rk45(X_model, x0, t0, ts, params);
  delta_x_1[1] =  x[1, 8] - x0[8] + 1e-5;
  for (i in 1:n_obs-1) {
    delta_x_1[i + 1] = x[i + 1, 8] - x[i, 8] + 1e-5;
  }
}
model {
  zeta ~ lognormal(0, 1);
  y1 ~ neg_binomial_2(delta_x_1, 5.54);
}
generated quantities {
  real log_lik;
  array[n_obs] int sim_y1;
  log_lik = neg_binomial_2_lpmf(y1 | delta_x_1, 5.54);
  sim_y1 = neg_binomial_2_rng(delta_x_1, 5.54);
}

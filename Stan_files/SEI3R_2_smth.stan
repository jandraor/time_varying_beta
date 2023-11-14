functions {
  vector X_model(real time, vector y, array[] real params) {
    vector[9] dydt;
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
    real adjust_Z_2;
    E_to_P = 0.73*0.3333333333*y[4];
    E_to_A = (1-0.73)*0.3333333333*y[4];
    P_to_I = 0.4761904762*y[5];
    A_to_R = 0.2*y[8];
    I_to_R = 0.3448275862*y[6];
    lambda = (params[1]*y[1]/4937796)*(y[6]+y[5]+0.5*y[8]);
    inv_nu = 1/params[3];
    P_to_C = P_to_I;
    S_to_E = lambda*y[3];
    adjust_Z = (y[2]-y[1])/(inv_nu/2);
    adjust_Z_2 = (params[2]-y[2])/(inv_nu/2);
    dydt[1] = adjust_Z;
    dydt[2] = adjust_Z_2;
    dydt[3] = -S_to_E;
    dydt[4] = S_to_E-E_to_P-E_to_A;
    dydt[5] = E_to_P-P_to_I;
    dydt[6] = P_to_I-I_to_R;
    dydt[7] = A_to_R+I_to_R;
    dydt[8] = E_to_A-A_to_R;
    dydt[9] = P_to_C;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  array[n_obs] int y1;
  real t0;
  array[n_obs] real ts;
}
parameters {
  real<lower = 0> zeta;
  real<lower = 0, upper = 1> upsilon;
  real<lower = 0, upper = 1> nu;
  real<lower = 0> P0;
}
transformed parameters{
  array[n_obs] vector[9] x; // Output from the ODE solver
  array[3] real params;
  vector[9] x0; // init values
  array[n_obs] real delta_x_1;
  x0[1] = 1; // Z
  x0[2] = 1; // Z_2
  x0[3] = (4937796) - P0; // S
  x0[4] = 0; // E
  x0[5] = P0; // P
  x0[6] = 0; // I
  x0[7] = 0; // R
  x0[8] = 0; // A
  x0[9] = 0; // C
  params[1] = zeta;
  params[2] = upsilon;
  params[3] = nu;
  x = ode_rk45(X_model, x0, t0, ts, params);
  delta_x_1[1] =  x[1, 9] - x0[9] + 1e-5;
  for (i in 1:n_obs-1) {
    delta_x_1[i + 1] = x[i + 1, 9] - x[i, 9] + 1e-5;
  }
}
model {
  zeta ~ lognormal(0, 1);
  upsilon ~ beta(2, 2);
  nu ~ beta(2, 2);
  P0 ~ lognormal(0, 1);
  y1 ~ poisson(delta_x_1);
}
generated quantities {
  real log_lik;
  array[n_obs] int sim_y1;
  log_lik = poisson_lpmf(y1 | delta_x_1);
  sim_y1 = poisson_rng(delta_x_1);
}

functions {
  vector SEI3R(real time, vector y, real[] params) {
    vector[15] dydt;
    real E_to_P;
    real P_to_A;
    real P_to_I;
    real A_to_R;
    real I_to_R;
    real lambda;
    real adjust_Z;
    real adjust_Z_2;
    real adjust_Z_3;
    real adjust_Z_4;
    real adjust_Z_5;
    real adjust_Z_6;
    real adjust_Z_7;
    real adjust_Z_8;
    real P_to_C;
    real S_to_E;
    E_to_P = 0.3333333333*y[10];
    P_to_A = (1-0.73)*0.4761904762*y[11];
    P_to_I = 0.73*0.4761904762*y[11];
    A_to_R = 0.2*y[14];
    I_to_R = 0.3448275862*y[12];
    lambda = (params[1]*y[1]/4937796)*(y[12]+y[11]+0.5*y[14]);
    adjust_Z = (y[2]-y[1])/(params[2]/8);
    adjust_Z_2 = (y[3]-y[2])/(params[2]/8);
    adjust_Z_3 = (y[4]-y[3])/(params[2]/8);
    adjust_Z_4 = (y[5]-y[4])/(params[2]/8);
    adjust_Z_5 = (y[6]-y[5])/(params[2]/8);
    adjust_Z_6 = (y[7]-y[6])/(params[2]/8);
    adjust_Z_7 = (y[8]-y[7])/(params[2]/8);
    adjust_Z_8 = (params[3]-y[8])/(params[2]/8);
    P_to_C = P_to_I;
    S_to_E = lambda*y[9];
    dydt[1] = adjust_Z;
    dydt[2] = adjust_Z_2;
    dydt[3] = adjust_Z_3;
    dydt[4] = adjust_Z_4;
    dydt[5] = adjust_Z_5;
    dydt[6] = adjust_Z_6;
    dydt[7] = adjust_Z_7;
    dydt[8] = adjust_Z_8;
    dydt[9] = -S_to_E;
    dydt[10] = S_to_E-E_to_P;
    dydt[11] = E_to_P-P_to_A-P_to_I;
    dydt[12] = P_to_I-I_to_R;
    dydt[13] = A_to_R+I_to_R;
    dydt[14] = P_to_A-A_to_R;
    dydt[15] = P_to_C;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  int<lower = 1> n_params;
  int<lower = 1> n_difeq;
  int y1[n_obs];
  real t0;
  real ts[n_obs];
}
parameters {
  real<lower = 0>            zeta;
  real<lower = 0, upper = 1> nu;
  real<lower = 0, upper = 1> upsilon;
  real<lower = 0>            P_0;
}
transformed parameters{
  vector[n_difeq] o[n_obs]; // Output from the ODE solver
  real y1_hat[n_obs];
  vector[n_difeq] y0;
  real params[n_params];
  y0[1] = 1;
  y0[2] = 1;
  y0[3] = 1;
  y0[4] = 1;
  y0[5] = 1;
  y0[6] = 1;
  y0[7] = 1;
  y0[8] = 1;
  y0[9] = 4937794 - P_0;
  y0[10] = 0;
  y0[11] = P_0;
  y0[12] = 0;
  y0[13] = 0;
  y0[14] = 0;
  y0[15] = 0;
  params[1] = zeta;
  params[2] = 1 / nu;
  params[3] = upsilon;
  o = ode_rk45(SEI3R, y0, t0, ts, params);
  y1_hat[1] =  o[1, 15]  - y0[15];
  for (i in 1:n_obs-1) {
    y1_hat[i + 1] = o[i + 1, 15] - o[i, 15] + 1e-4;
  }
}
model {
  zeta    ~ lognormal(0, 1);
  upsilon ~ beta(2, 2);
  nu      ~ beta(2, 2);
  P_0     ~ lognormal(0, 1);
  y1      ~ poisson(y1_hat);
}
generated quantities {
  real log_lik;
  log_lik = poisson_lpmf(y1 | y1_hat);
}

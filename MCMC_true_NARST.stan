// I think the function should go here, perhaps? 
functions{
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
    real H = y[1];
    real A = y[2];
    real E = y[3];
    // real N = x_i[1]; I hope I don't need this? Because...I don't have N?
    
    real gamma = theta[1]; 
    real LAMBDA_H = theta[2]; 
    real beta_HH = theta[3]; 
    real beta_AH = theta[4]; 
    real beta_EH = theta[5];
    real LAMBDA_A = theta[6];
    real beta_HA = theta[7]; 
    real beta_AA = theta[8]; 
    real beta_EA = theta[9]; 
    real LAMBDA_E = theta[10]; 
    real beta_HE = theta[11]; 
    real beta_AE = theta[12]; 
    real beta_EE = theta[13]; 
    real mu = theta[14];
    
    
    real dH_dt = gamma * LAMBDA_H * (1-H) + LAMBDA_H * beta_HH * H * (1-H) + LAMBDA_H * beta_AH * (1-H) * A + LAMBDA_H * beta_EH * (1-H) * E - mu * H;
    real dA_dt = gamma * LAMBDA_A * (1-A) + LAMBDA_A * beta_AA * A * (1-A) + LAMBDA_A * beta_HA * (1-A) * H + LAMBDA_A * beta_EA * (1-A) * E - mu * A;
    real dE_dt = LAMBDA_E * beta_EE * E * (1-E) + LAMBDA_E * beta_HE * (1-E) * H + LAMBDA_E * beta_AE * (1-E) * A - mu * E;
    
    return{dH_dt, dA_dt, dE_dt};
  }
} 


// Fixed data is declared in the data block
data {
  int<lower=1> n_days;
  real y0[3]; 
  real t0; 
  real ts[n_days];
  int human_prev[n_days];
}

// Then we transform the data, in this case to match our signature of sir
// This block is evaluated once
transformed data {
  real x_r[0];
  int x_i[0]; //of length 0 because we have nothing to put in it (?)
}


// The parameters accepted by the model
parameters {
  real<lower=0, upper=1> gamma; 
  real<lower=0, upper=1> LAMBDA_H; 
  real<lower=0, upper=1> beta_HH;
  real<lower=0, upper=1> beta_AH;
  real<lower=0, upper=1> beta_EH; 
  real<lower=0, upper=1> LAMBDA_A;
  real<lower=0, upper=1> beta_HA; 
  real<lower=0, upper=1> beta_AA; 
  real<lower=0, upper=1> beta_EA; 
  real<lower=0, upper=1> LAMBDA_E; 
  real<lower=0, upper=1> beta_HE; 
  real<lower=0, upper=1> beta_AE; 
  real<lower=0, upper=1> beta_EE; 
  real<lower=0, upper=1> mu;
  real<lower=0> phi_inv;
}

// And then transform the parameters...
// This block is evaluated/differentiated at each leapfrog step (multiple times per iteration)
// Main computational cost! Sadly... 
transformed parameters{
  real y[n_days, 3];
  real phi = 1. / phi_inv;
  {
    real theta[14];
    theta[1] = gamma; 
    theta[2] = LAMBDA_H; 
    theta[3] = beta_HH; 
    theta[4] = beta_AH; 
    theta[5] = beta_EH; 
    theta[6] = LAMBDA_A; 
    theta[7]= beta_HA; 
    theta[8] = beta_AA; 
    theta[9] = beta_EA; 
    theta[10] = LAMBDA_E; 
    theta[11] = beta_HE; 
    theta[12] = beta_AE;  
    theta[13] = beta_EE; 
    theta[14] = mu; 
  
    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}

// The model to be estimated
model {
  // First, the priors
  // All of these are bounded between 0 and 1! 
  // Basic assumption: mean = 0.5, sd = 0.1
  // I know normal is not correct -- runif seems like an option, but not sure how to input restrictions/assumptions with runif
  gamma ~ normal(0.5, 0.1); 
  LAMBDA_H ~ normal(0.5, 0.1); 
  beta_HH ~ normal(0.6, 0.1); // Assumption: beta_HH > beta_AH, beta_EH
  beta_AH ~ normal(0.4, 0.1); // Ditto above assumption
  beta_EH ~ normal(0.3, 0.1); // Assumption: beta_EA > beta_EH
  LAMBDA_A ~ normal(0.7, 0.1); // Assumption: lambda_A > lambda_H, majority global ABU in animals raised for food (73%) (Booton et al.)
  beta_HA ~ normal(0.5, 0.1);
  beta_AA ~ normal(0.6, 0.1); // Assumption: beta_AA > beta_HA, beta_EA
  beta_EA ~ normal(0.5, 0.1);
  LAMBDA_E ~ normal(0.4, 0.1); // Assumption: lambda_H > lambda_E
  beta_HE ~ normal(0.5, 0.1);
  beta_AE ~ normal(0.5, 0.1);
  beta_EE ~ normal(0.5, 0.1);
  mu ~ normal(0.5, 0.1);
  
  // Second, the sampling distributions
  // col(matrix x, int n) - the n-th column of matrix x 
  human_prev ~ neg_binomial_2(col(to_matrix(y), 1), phi);
  // where col(to_matrix(y), 1) is the first column of ODE solutions - human infections
}

// We can make a block of generated quantities
// This block is evaluated once per iteration
generated quantities {
  real pred_human_prev[n_days];
  pred_human_prev = neg_binomial_2_rng(col(to_matrix(y), 1), phi);
  
}




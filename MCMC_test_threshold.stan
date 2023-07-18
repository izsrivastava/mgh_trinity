//
// This Stan program defines a simple MCMC model 
//


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


// The input data is a vector 'y' of length 'N' (?)
// Fixed data is declared in the data block
data {
  int<lower=1> n_days;
  real y0[3]; 
  real t0; 
  real ts[n_days];
  int human_prev[n_days];
  int animal_prev[n_days];
  int enviro_prev[n_days];
}

// Then we transform the data, in this case to match our signature of sir
// This block is evaluated once
transformed data {
  real x_r[0];
  int x_i[0]; //of length 0 because we have nothing to put in it (?)
}


// The parameters accepted by the model
parameters {
  real<lower=0.3, upper=0.751> gamma; 
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
  gamma ~ uniform(0.3, 0.751); // Range: 30-75.1%
  LAMBDA_H ~ beta(4,2); // Estimate: high avg. in 3 months=37.7%, so want most data below that
  beta_HH ~ beta(2,2); // Assumption: >BETA_AH // Estimate from Booton: 50%
  beta_AH ~ uniform(0,1); // No estimate // Assumption: >BETA_HA, <BETA_HH, <BETA_AA
  beta_EH ~ beta(2,4); // Assumption: <BETA_EA, <BETA_HE // Estimate: 28.5-31% 
  LAMBDA_A ~ beta(5.5, 2); // Estimate from literature review
  beta_HA ~ uniform(0,1); // No estimate // Assumption: <BETA_AH, <BETA_AA, <BETA_HH 
  beta_AA ~ beta(4,2); // Estimate from Booton // Assumption: >BETA_HA, >BETA_AH
  beta_EA ~ beta(2,1); // Assumption: >BETA_EH // Estimate: 65% 
  LAMBDA_E ~ beta(3,2); // Range: 0.45 * 30-90% of LAMBDA_A (estimated that at 60%)
  beta_HE ~ beta(3,2); // Assumption: >BETA_EH // Estimate: 66% 
  beta_AE ~ beta(3,2); // Estimate: 66% 
  beta_EE ~ beta(4,2); // Assumption: >BETA_HE, >BETA_EH // Estimate: 48-84%
  mu ~ uniform(0,1); // No assumptions/estimates, full parameter space
  
  // Second, the sampling distributions
  // col(matrix x, int n) - the n-th column of matrix x 
  human_prev ~ neg_binomial_2(col(to_matrix(y), 1), phi);
  animal_prev ~ neg_binomial_2(col(to_matrix(y), 2), phi);
  enviro_prev ~ neg_binomial_2(col(to_matrix(y), 3), phi);
  // where col(to_matrix(y), 1) is the first column of ODE solutions - human infections
}

// We can make a block of generated quantities
// This block is evaluated once per iteration
generated quantities {
  real pred_human_prev[n_days];
  pred_human_prev = neg_binomial_2_rng(col(to_matrix(y), 1), phi);
}

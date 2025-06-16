functions {
  vector compute_p_k(real p, int K) {
    vector[K] p_k;
    for (k in 1:K)
      p_k[k] = pow(1 - p, k - 1);
    return p_k;
  }
}

data {
  int<lower=1> I;                 // number of age groups
  int<lower=1> J;                 // number of weeks
  int<lower=0> Y[I, J];           // observed reported cases
  vector[J] ratio;                // R_j: external reproduction number
  vector<lower=0>[I] mean_X0;        // prior mean for initial latent infections
  vector<lower=0>[I] sd_X0;       // prior SD for initial infections
  real<lower=0> phi_p;            // exponential rate for pi random walk
  
  // Hospitalization model
  int<lower=1> K;                     // convolution window size
  int<lower=0> H[I, J];              // observed hospitalizations
  real<lower=0> alpha;              // shape for geometric delay
  real<lower=0> beta;               // rate for geometric delay
  real phi_tau;   
}

parameters {
  matrix<lower=0>[I, J] X;           // latent infections
  vector[J] logit_pi;                // logit of shared reporting probability pi_j
  real<lower=0> psi_p;               // log sd of pi RW(1)
  // matrix<lower=0>[I, I] C_raw;    // raw (symmetric) contact matrix
  vector<lower=0>[(I*(I+1)/2)] C_raw; // I * (I + 1) / 2
  
  // Hospitalization model
  matrix<lower=0>[I, J] A;           // latent hospital admissions
  real<lower=0, upper=1> theta;      // geometric delay param
  matrix[I, J] tau_logit;            // admission age-time probabilities. 
  real<lower=0> psi_tau;             // AR(1) scale
  real<lower=-1, upper=1> rho;       // Autoregressive correlation parameter
}

transformed parameters {
  vector<lower=0, upper=1>[J] pi;
  real sd_p = exp(-0.5 * psi_p);
  matrix[I, I] C;


  for (j in 1:J)
    pi[j] = inv_logit(logit_pi[j]);

    {
    int idx = 1;
    for (i in 1:I) {
      for (j in i:I) {
        C[i, j] = C_raw[idx];
        C[j, i] = C_raw[idx];
        idx += 1;
      }
    }
  }
  
  // Hospital model
  vector[K] p_k = compute_p_k(theta, K);            // delay probabilities
  matrix<lower=0, upper=1>[I, J] tau;           // admission probabilities
  matrix<lower=0>[I, J] H_mean;                 // Hospital data mean
  real<lower=0> sd_tau = exp(-0.5*psi_tau);

  for (i in 1:I)
    for (j in 1:J)
      tau[i, j] = inv_logit(tau_logit[i, j]);
  
   for (i in 1:I) {
    for (j in 1:J) {
      real mu = 0;
      for (k in 1:min(K, j)) {
        mu += A[i, j - k + 1] * p_k[k];
      }
      H_mean[i, j] = mu;
    }
  }
    
}

model {
  
// Infection model
  
  // Prior on logit(pi)
  logit_pi[1] ~ normal(0, 1);
  for (j in 2:J)
    logit_pi[j] ~ normal(logit_pi[j - 1], sd_p);
    
  // Prior on initial latent infections
  for (i in 1:I)
    X[i, 1] ~ normal(mean_X0[i], sd_X0[i]);
    
    // Infection dynamics
  for (j in 2:J) {
    for (i in 1:I) {
      real lambda_ij = 0;
      for (ip in 1:I)
        lambda_ij += C[i, ip] * X[ip, j - 1];
      lambda_ij *= ratio[j];
      X[i, j] ~ normal(lambda_ij, sqrt(lambda_ij));
    }
  }
  
    // Observation model
  for (i in 1:I)
    for (j in 1:J)
      Y[i, j] ~ normal(X[i, j]*pi[j], sqrt(X[i, j]*pi[j]*(1-pi[j])));
  
  
  // priors
  // Prior on upper triangle
  to_vector(C_raw) ~ normal(1.5,1);  
      
   target += log(0.5 * phi_p) - phi_p*exp(-0.5*psi_p) - 0.5*psi_p;
   
// hospitalization process
   
  // Prior on logit(tau)
  for (i in 1:I) {
  tau_logit[i, 1] ~ normal(0,1);
  for (j in 2:J) {
    tau_logit[i, j] ~ normal(rho * tau_logit[i, j - 1], sd_tau);
  }
}

  target += log(0.5 * phi_tau) - phi_tau*exp(-0.5*psi_tau) - 0.5*psi_tau;
  target += gamma_lpdf(theta | alpha, beta);         //prior for geometric delay parameter
  
  for (i in 1:I){
    for (j in 1:J){
      target += normal_lpdf(A[i,j]| X[i,j]*tau[i,j], sqrt(X[i,j]*tau[i,j]*(1-tau[i,j])));
    }
  }
  
  for (i in 1:I){
    for (j in 1:J){
      target+=  poisson_lpmf(H[i, j]|H_mean[i, j]);}}

  // Priors for AR(1) hyperparameters
  rho ~ uniform(-1, 1);

}

generated quantities {
  vector[J] X_total;
  for (j in 1:J)
    X_total[j] = sum(X[, j]);
}





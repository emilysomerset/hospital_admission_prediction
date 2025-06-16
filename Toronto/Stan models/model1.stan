functions {
  vector compute_p_k(real p, int K) {
    vector[K] p_k;
    for (k in 1:K)
      p_k[k] = pow(1 - p, k - 1);
    return p_k;
  }

  vector logit_inverse_rw(vector x) {
    int J = num_elements(x);
    vector[J] out;
    for (j in 1:J)
      out[j] = inv_logit(x[j]);
    return out;
  }
}


data {
  int<lower=1> I;                     // number of age groups
  int<lower=1> J;                     // number of weeks
  int<lower=1> K;                     // convolution window size
  
  // Infection Model
  int<lower=0> Y[I, J];              // observed reported cases by age
  int<lower=0> case_counts[J];   //observed reported cases (the sum of each column of Y)
  vector[J] ratio;                  // wastewater-reproduction number
  real mean_z0;                     // prior mean for initial infection count
  real sd_z0;                       // prior SD for initial infection count
  real phi_p;   
  
  // Hospital model
  int<lower=0> H[I, J];              // observed hospitalizations
  real<lower=0> alpha;              // shape for geometric delay
  real<lower=0> beta;               // rate for geometric delay
  real phi_tau;   
}

parameters {
  // Infection model
  vector[J] Z_raw;                   // latent infection counts
  vector[J] logit_p_std;             // standardized random walk steps
  real logit_p0;                     // initial logit(p) value
  real psi_p;                        // RW(1) scale parameter

  
  // Hospitalization model
  matrix<lower=0>[I, J] A;           // latent hospital admissions
  real<lower=0, upper=1> theta;             // geometric delay param
  matrix[I, J] tau_std;  // standard normal innovations
  vector[I] tau_logit0;
  real<lower=0> psi_tau;             // AR(1) scale
  real<lower=-1, upper=1> rho;         // Autoregressive correlation parameter
}

transformed parameters {

  // Infection model
  vector[J] p;
  vector[J] pZ;
  vector[J] ppZ;
  real<lower=70> Z[J];
  real logit_p[J];  // reporting probability
  real<lower=0> sd_p = exp(-0.5*psi_p);
  
  logit_p[1] = logit_p0;
  for (j in 2:J){
  logit_p[j] = logit_p[j - 1] + logit_p_std[j] * sd_p;}
  
  Z[1] = mean_z0 + sd_z0*Z_raw[1];
  for (i in 2:J){
    Z[i] = ratio[i]*Z[i-1] + sqrt(ratio[i]*Z[i-1])*Z_raw[i];
  }
  
  for (i in 1:J){
    p[i] = exp(logit_p[i])/(1+exp(logit_p[i]));
    pZ[i] = p[i]*Z[i]; //pZ(0) = pZ1 = p1*Z1
    ppZ[i] = pow(p[i]*(1-p[i])*Z[i],0.5);
  }
  
  // Hospital model
  vector[K] p_k = compute_p_k(theta, K);            // delay probabilities
  matrix[I, J] tau_logit;
  matrix<lower=0, upper=1>[I, J] tau;           // admission probabilities
  matrix<lower=0>[I, J] H_mean;                 // Hospital data mean
  real<lower=0> sd_tau = exp(-0.5*psi_tau);
  
  for (i in 1:I) {
  tau_logit[i, 1] = tau_logit0[i];
  for (j in 2:J) {
    tau_logit[i, j] = rho * tau_logit[i, j - 1] + tau_std[i, j] * sd_tau;
  }
}

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
  for (i in 1:J){
   target += normal_lpdf(case_counts[i]|pZ[i], sqrt(ppZ[i]));}
   target += normal_lpdf(logit_p0|0,1);
   target += normal_lpdf(logit_p_std|0,1);
   target += normal_lpdf(Z_raw|0,1);
   target += log(0.5 * phi_p) - phi_p*exp(-0.5*psi_p) - 0.5*psi_p;
  
  // Hospital model
  target += log(0.5 * phi_tau) - phi_tau*exp(-0.5*psi_tau) - 0.5*psi_tau;
  target += gamma_lpdf(theta | alpha, beta);         //prior for geometric delay parameter
  
  for (i in 1:I){
    for (j in 1:J){
      target += normal_lpdf(A[i,j]| Y[i,j]*tau[i,j], sqrt(Y[i,j]*tau[i,j]+1e-2));
    }
  }
  
  for (i in 1:I){
    for (j in 1:J){
      target+=  poisson_lpmf(H[i, j]|H_mean[i, j]);}}
  
  target += normal_lpdf(to_vector(tau_std) | 0, 1);
  target += normal_lpdf(tau_logit0 | 0, 1);


  // Priors for AR(1) hyperparameters
  rho ~ uniform(-1, 1);
}





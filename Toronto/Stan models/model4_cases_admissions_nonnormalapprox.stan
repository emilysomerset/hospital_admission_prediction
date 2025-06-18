functions {
real multinomial_log_approx(int A, int Y, real X, vector p) {
      // validity check first
    if (X < Y || A < 0 || A > Y) {
      return negative_infinity();
    }
    
  real x1 = A;
  real x2 = Y - A;
  real x3 = X - Y;
  real log_pmf = lgamma(X + 1)
                 - lgamma(x1 + 1)
                 - lgamma(x2 + 1)
                 - lgamma(x3 + 1)
                 + x1 * log(p[1])
                 + x2 * log(p[2])
                 + x3 * log(p[3]);
    return log_pmf;
}
}

data {
  int<lower=1> I;                 // number of age groups
  int<lower=1> J;                 // number of weeks
  int<lower=0> Y[I, J];           // observed reported cases
  vector[J] ratio;                // R_j: external reproduction number
  vector<lower=0>[I] mean_X0;     // prior mean for initial latent infections
  vector<lower=0>[I] sd_X0;       // prior SD for initial infections
  real<lower=0> phi_p;            // exponential rate for pi random walk
  
  // Admissions
  int<lower=0> A[I,J];            // observed hospital admissions for COVID-19
  real phi_tau;                   // exponential rate for tau autoregressive
}

parameters {
  matrix<lower=0>[I, J] X;        // latent infections
  real logit_pi_init;
  vector[J - 1] z_logit_pi;      // standardized logit of shared reporting probability pi_j
  real<lower=0> psi_p;            // log sd of pi RW(1)
  // matrix<lower=0>[I, I] C_raw;    // raw (symmetric) contact matrix
  vector<lower=0>[(I*(I+1)/2)] C_raw; // I * (I + 1) / 2
  
  // Admissions
  matrix[I, J] tau_logit;            // admission age-time probabilities. 
  real<lower=0> psi_tau;             // AR(1) scale
  // real<lower=-1, upper=1> rho;       // Autoregressive correlation parameter
}

transformed parameters {
  real sd_p = exp(-0.5 * psi_p);
  vector[J] logit_pi;
  
  logit_pi[1] = logit_pi_init;
  for (j in 2:J)
    logit_pi[j] = logit_pi[j - 1] + sd_p * z_logit_pi[j - 1];
    
  matrix[I, I] C;
    
    
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
  
  // Admissions
  real<lower=0> sd_tau = exp(-0.5*psi_tau);
  array[I, J] simplex[3] p;
  matrix[I,J] denom;
  
  for (i in 1:I){
  for (j in 1:J){
    denom[i,j] = 1 + exp(tau_logit[i,j]) + exp(logit_pi[j]);
    p[i, j][1] = exp(tau_logit[i, j]) / denom[i,j];
    p[i, j][2] = exp(logit_pi[j]) / denom[i,j];
    p[i, j][3] = 1 / denom[i,j];
  }
  }
    
}

model {
  
  // Prior on logit(pi)
  logit_pi_init ~ normal(0, 1);
  z_logit_pi ~ normal(0, 1);
    
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
  
   // Prior on logit(tau)
  for (i in 1:I) {
  tau_logit[i, 1] ~ normal(-2,1);
  for (j in 2:J) {
    tau_logit[i, j] ~ normal(tau_logit[i, j - 1], sd_tau);
  }
}

for (i in 1:I) {
    for (j in 1:J) {
      target += multinomial_log_approx(A[i, j], Y[i, j], X[i, j], p[i,j]);
    }
  }
  
  
  // priors
  // Prior on upper triangle
  to_vector(C_raw) ~ normal(1.5,1);  
      
   target += log(0.5 * phi_p) - phi_p*exp(-0.5*psi_p) - 0.5*psi_p;
  
  //admissions 
  target += log(0.5 * phi_tau) - phi_tau*exp(-0.5*psi_tau) - 0.5*psi_tau;
  // rho ~ uniform(-1, 1);         //prior for geometric delay parameter

}

generated quantities {
  vector[J] X_total;
  for (j in 1:J)
    X_total[j] = sum(X[, j]);
}
data {
  int<lower=0> N;                  //number of data points
  int<lower=0> lag;
  real<lower=0> case_counts[N-lag];               //weekly case counts
  vector[N] ratio;     //geometric derivative wastewater
  int<lower=0> obs_start_case;
  real phi1;         
  real mean_z0;
  real sd_z0;
}

parameters {
  vector[N-lag] Z_raw;
  real logit_p_raw[N-lag];  // reporting probability
  real theta_p;
}

transformed parameters {
  vector[N-lag] p;
  vector[N-lag] pZ;
  vector[N-lag] ppZ;
  real<lower=0> Z[N-lag];
  real logit_p[N-lag];  // reporting probability
  real<lower=0> sd_p;
  
  sd_p = exp(-0.5*theta_p);
  
  logit_p[1] = logit_p_raw[1];
  for (i in 2:(N-lag)){
    logit_p[i] = logit_p[i-1] + logit_p_raw[i]*sd_p;
  }
  
  Z[1] = mean_z0 + sd_z0*Z_raw[1];
  for (i in 2:(N-lag)){
    Z[i] = ratio[i]*Z[i-1] + sqrt(ratio[i]*Z[i-1])*Z_raw[i];
  }
  
  for (i in 1:(N-lag)){
    p[i] = exp(logit_p[i])/(1+exp(logit_p[i]));
    pZ[i] = p[i]*Z[i]; //pZ(0) = pZ1 = p1*Z1
    ppZ[i] = pow(p[i]*(1-p[i])*Z[i],0.5);
  }
  
}

model {
  for (i in 1:(N-lag)){
   target += normal_lpdf(case_counts[i+lag]|pZ[i], sqrt(ppZ[i]));}
   
   //priors
  target += normal_lpdf(logit_p_raw|0,1);
  target += normal_lpdf(Z_raw|0,1);
  target += log(0.5 * phi1) - phi1*exp(-0.5*theta_p) - 0.5*theta_p;
}





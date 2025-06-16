#include <TMB.hpp>
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  
  // Data
  DATA_INTEGER(I);
  DATA_INTEGER(J);
  DATA_MATRIX(Y);         // [I x J]
  DATA_MATRIX(A);         // [I x J]
  DATA_VECTOR(ratio);     // [J]
  DATA_VECTOR(mean_X0);   // [I]
  DATA_VECTOR(sd_X0);     // [I]
  DATA_SCALAR(phi_p);
  DATA_SCALAR(phi_tau);
  
  // Parameters
  PARAMETER_MATRIX(X);               // [I x J]
  PARAMETER(logit_pi_init);
  PARAMETER_VECTOR(z_logit_pi);      // [J - 1]
  PARAMETER(psi_p);
  PARAMETER_VECTOR(C_raw);           // length = I*(I+1)/2
  PARAMETER_MATRIX(tau_logit);       // [I x J]
  PARAMETER(psi_tau);
  PARAMETER(rho);                    // AR(1)
  
  // Transformed parameters
  Type sd_p = exp(-0.5 * psi_p);
  vector<Type> logit_pi(J);
  logit_pi(0) = logit_pi_init;
  for (int j = 1; j < J; ++j) {
    logit_pi(j) = logit_pi(j - 1) + sd_p * z_logit_pi(j - 1);
  }
  
  // Build symmetric contact matrix C
  matrix<Type> C(I, I);
  int idx = 0;
  for (int i = 0; i < I; ++i) {
    for (int j = i; j < I; ++j) {
      C(i,j) = C_raw(idx);
      C(j,i) = C_raw(idx);
      idx++;
    }
  }
  
  Type nll = 0.0;
  
  // Prior on logit_pi
  nll -= dnorm(logit_pi_init, Type(0.0), Type(1.0), true);
  for (int j = 0; j < J - 1; ++j)
    nll -= dnorm(z_logit_pi(j), Type(0.0), Type(1.0), true);
  
  // Prior on initial X
  for (int i = 0; i < I; ++i) {
    nll -= dnorm(X(i, 0), mean_X0(i), sd_X0(i), true);
  }
  
  // Infection dynamics
  for (int j = 1; j < J; ++j) {
    for (int i = 0; i < I; ++i) {
      Type lambda = 0;
      for (int ip = 0; ip < I; ++ip)
        lambda += C(i, ip) * X(ip, j - 1);
      lambda *= ratio(j);
      nll -= dnorm(X(i, j), lambda, sqrt(lambda), true);
    }
  }
  
  // AR(1) prior for tau_logit
  Type sd_tau = exp(-0.5 * psi_tau);
  for (int i = 0; i < I; ++i) {
    nll -= dnorm(tau_logit(i,0), Type(-2), Type(1), true);
    for (int j = 1; j < J; ++j) {
      Type mean = rho * tau_logit(i,j - 1);
      nll -= dnorm(tau_logit(i,j), mean, sd_tau, true);
    }
  }
  
  // Likelihood (multivariate normal approximation to multinomial)
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j < J; ++j) {
      Type exp_tau = exp(tau_logit(i,j));
      Type exp_pi = exp(logit_pi(j));
      Type denom = 1 + exp_tau + exp_pi;
      vector<Type> p(3);
      p(0) = exp_tau / denom;
      p(1) = exp_pi / denom;
      p(2) = 1.0 / denom;
      
      Type total_X = X(i,j);
      vector<Type> mu = total_X * p;
      
      vector<Type> y_vec(3);
      y_vec(0) = A(i,j);
      y_vec(1) = Y(i,j) - A(i,j);
      y_vec(2) = X(i,j) - Y(i,j);
      
      matrix<Type> Sigma(3,3);
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          Sigma(k,l) = -total_X * p(k) * p(l);
      for (int k = 0; k < 3; ++k)
        Sigma(k,k) += total_X * p(k);
      
      for (int k = 0; k < 3; ++k)
        Sigma(k,k) += 1e-6;  // stabilization
      
      nll += MVNORM(Sigma)(y_vec - mu);
    }
  }
  
  // Prior on contact matrix elements
  for (int i = 0; i < C_raw.size(); ++i)
    nll -= dnorm(C_raw(i), Type(1.5), Type(1.0), true);
  
  // Half-Cauchy(0, phi_p) prior on sd_p (implemented via psi_p)
  nll -= log(0.5 * phi_p) - phi_p * sd_p - 0.5 * psi_p;
  nll -= log(0.5 * phi_tau) - phi_tau * sd_tau - 0.5 * psi_tau;
  
  // rho prior: uniform(-1,1)
  if (rho < -1.0 || rho > 1.0) nll += INFINITY;
  
  return nll;
}

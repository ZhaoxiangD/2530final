//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of years
  int<lower=0> K; // number of age group
  int<lower=0> M; // number of cancer types
  matrix[M, K] y13; // cancer incidence for 2013
  matrix[M, K] y14; // cancer incidence for 2014
  matrix[M, K] y15; // cancer incidence for 2015
  matrix[M, K] y16; // cancer incidence for 2016
  matrix[M, K] y17; // cancer incidence for 2017
  matrix[K, N] pop; // population
  matrix[K, K] X; // age group indicator matrix
  matrix[M, M] C; // cancer type corrlation matrix
  real tar_sigma1;
  real tar_sigma2;
  real rho1;
  real rho2;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower = 0.0001> inv_logit_p[M,N,K];
  real<lower = 0.0001> log_theta[M,N,K];
  real <lower = 0> epi1;
  real <lower = 0> epi2;
  vector [K] alpha1;
  vector [K] alpha2;
  vector [K] beta;
  vector [M] gamma1;
  vector [M] gamma2;
  real sigma1_initial;
  real sigma2_initial;
  real epsilon1[M, N];
  real epsilon2[M, N];
  real k1;
  real k2;
  real<lower = 0.00001> inv_tau_gamma1;
  real<lower = 0.0001> inv_tau_gamma2;
  real<lower = 0.0001> inv_tau_epsilon2[M,N];
}

transformed parameters {
  real<lower = 0> p[M,N,K] = inv_logit(inv_logit_p);
  real theta[M,N,K] = exp(log_theta);
  real<lower = 0> tau_gamma1 = 1/(inv_tau_gamma1)^2;
  real<lower = 0> tau_gamma2 = 1/(inv_tau_gamma2)^2;
  real tau_epsilon2[M,N];
  vector [M] mu;
  real sigma1[N];
  real sigma2[N];
  
  
  for (m in 1:M){
    mu[m] = 0;
  }
  
  for (m in 1:M){
    for (n in 1:N){
      tau_epsilon2[m,n] = 1/(inv_tau_epsilon2[m,n])^2;
    }
  }
  
  for (n in 1:N){
    if (n ==1){
      sigma1[n] = sigma1_initial;
      sigma2[n] = sigma2_initial;
    } else {
      sigma1[n] = sigma1[n-1] * rho1 + k1;
      sigma2[n] = sigma2[n-1] * rho2 + k2;
    }
}
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model{
  alpha1 ~ normal(0, 0.001);
  alpha2 ~ normal(0, 0.001);
  beta ~ normal(0, 0.001);
  inv_tau_gamma1 ~ cauchy(0, 100);
  inv_tau_gamma2 ~ cauchy(0, 100);
  gamma1 ~ multi_normal(mu, tau_gamma1*C);
  gamma2 ~ multi_normal(mu, tau_gamma2*C);
  
  epi1 ~ cauchy(0, 100);
  epi2 ~ cauchy(0, 100);
  
  for (m in 1:M){
    for (n in 1:N){
      inv_tau_epsilon2[m,n] ~ cauchy(0, 100);
      epsilon2[m,n] ~ normal(0, inv_tau_epsilon2[m,n]);
    }
  }
  
  sigma1_initial ~ normal(0, tar_sigma1);
  sigma2_initial ~ normal(0, tar_sigma2);
  
  k1 ~ normal(0, tar_sigma1);
  k2 ~ normal(0, tar_sigma2);
  
  for (m in 1:M){
    for (k in 1:K){
      for (n in 1:N){
        inv_logit_p[m,n,k] ~ normal(col(X,k)' * alpha1 + log(pop[k,n])* col(X,k)' * alpha2 + gamma1[m] + sigma1[n], epi1);
        log_theta[m,n,k] ~ normal(col(X,k)' * beta + log(pop[k,n]) + gamma2[m] + sigma2[n] + epsilon2[m,n], epi2);
      }
    }
  }
  
  // for (k in 1:K){
  //   for (n in 1:N){
  //     vector [M] mu_p1;
  //     vector [M] mu_theta1;
  //     for (m in 1:M){
  //       mu_p1[m] = col(X,k)' * alpha1 + log(pop[k,n])* col(X,k)' * alpha2;
  //       mu_theta1[m] = col(X,k)' * beta + log(pop[k,n]);
  //     }
  //     to_vector(inv_logit_p[,n,k]) ~ multi_normal(mu_p1 + sigma2[n], tau_gamma1*C);
  //     to_vector(log_theta[,n,k]) ~ multi_normal(mu_theta1 + sigma2[n], tau_gamma2*C);
  //   }
  // }
  // 
  // 
  // for (k in 1:K){
  //   for (m in 1:M){
  //     for (n in 1:N){
  //       real  mu_theta2 = log_theta[m,n,k];
  //       log_theta[m,n,k] ~ normal(mu_theta2, tau_epsilon2[m,n]);
  //     }
  //   }
  // }
  // 
  // 2013
  for (m in 1:M){
    for (k in 1:K){
      if (y13[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,1,k]);
      } else {
        target += log(p[m,1,k]) + poisson_lpmf(to_int(y13[m,k]) | theta[m,1,k]);
      }
    }
  }
  
  // 2014
  for (m in 1:M){
    for (k in 1:K){
      if (y14[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,2,k]);
      } else {
        target += log(p[m,2,k]) + poisson_lpmf(to_int(y14[m,k]) | theta[m,2,k]);
      }
    }
  }
  
  // 2015
  for (m in 1:M){
    for (k in 1:K){
      if (y15[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,3,k]);
      } else {
        target += log(p[m,3,k]) + poisson_lpmf(to_int(y15[m,k]) | theta[m,3,k]);
      }
    }
  }

  // 2016
  for (m in 1:M){
    for (k in 1:K){
      if (y16[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,4,k]);
      } else {
        target += log(p[m,4,k]) + poisson_lpmf(to_int(y16[m,k]) | theta[m,4,k]);
      }
    }
  }
  
  // 2017
  for (m in 1:M){
    for (k in 1:K){
      if (y17[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,5,k]);
      } else {
        target += log(p[m,5,k]) + poisson_lpmf(to_int(y17[m,k]) | theta[m,5,k]);
      }
    }
  }
}


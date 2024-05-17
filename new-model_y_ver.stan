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
  real<lower = 0> tau_gamma1;
  real<lower = 0> tau_gamma2;
  real<lower = 0> tau_epsilon2;


}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
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
}

transformed parameters {
  vector [M] mu;
  real sigma1[N];
  real sigma2[N];
  
  for (m in 1:M){
    mu[m] = 0;
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
  alpha1 ~ normal(0, 1);
  alpha2 ~ normal(0, 1);
  beta ~ normal(0, 1);
  gamma1 ~ multi_normal(mu, tau_gamma1*C);
  gamma2 ~ multi_normal(mu, tau_gamma2*C);

  real inv_logit_p[M,N,K];
  real log_theta[M,N,K];
  
  for (m in 1:M){
    for (n in 1:N){
      epsilon2[m,n] ~ normal(0, tau_epsilon2);
    }
  }
  
  sigma1_initial ~ normal(0, tar_sigma1);
  sigma2_initial ~ normal(0, tar_sigma2);
  
  k1 ~ normal(0, tar_sigma1);
  k2 ~ normal(0, tar_sigma2);
  
  for (m in 1:M){
    for (k in 1:K){
      for (n in 1:N){
        inv_logit_p[m,n,k] = col(X,k)' * alpha1 + log(pop[k,n])* col(X,k)' * alpha2 + gamma1[m] + sigma1[n];
        log_theta[m,n,k]  = col(X,k)' * beta + log(pop[k,n]) + gamma2[m] + sigma2[n] + epsilon2[m,n];
      }
    }
  }

  real p[M,N,K] = inv_logit(inv_logit_p);
  real theta[M,N,K] = exp(log_theta);

  
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
  
  matrix[M,M] v_mat = tau_gamma1*C;
  
  // 2013
  for (m in 1:M){
    for (k in 1:K){
      if (y13[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,1,k]) + normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[1] | 0, tar_sigma1);
      } else {
        target += log(p[m,1,k]) + poisson_lpmf(to_int(y13[m,k]) | theta[m,1,k])+ normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[1] | 0, tar_sigma1) +normal_lpdf(gamma2[m] | 0,v_mat[m,m]) + normal_lpdf(sigma2[1] | 0, tar_sigma2) + normal_lpdf(beta[k] |0,1) + normal_lpdf(epsilon2[m,1] | 0, tau_epsilon2);
      }
    }
  }
  
  // 2014
  for (m in 1:M){
    for (k in 1:K){
      if (y14[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,2,k]) + normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[2] | sigma1[1]* rho1, tar_sigma1);
      } else {
        target += log(p[m,2,k]) + poisson_lpmf(to_int(y14[m,k]) | theta[m,2,k])+ normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[2] | sigma1[1]* rho1, tar_sigma1) +normal_lpdf(gamma2[m] | 0,v_mat[m,m]) + normal_lpdf(sigma2[2] |sigma2[1]* rho1, tar_sigma2) + normal_lpdf(beta[k] |0,1) + normal_lpdf(epsilon2[m,2] | 0, tau_epsilon2);
      }
    }
  }
  
  // 2015
  for (m in 1:M){
    for (k in 1:K){
      if (y15[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,3,k]) + normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[3] | sigma1[2]* rho1, tar_sigma1);
      } else {
        target += log(p[m,3,k]) + poisson_lpmf(to_int(y15[m,k]) | theta[m,3,k])+ normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[3] | sigma1[2]* rho1, tar_sigma1) +normal_lpdf(gamma2[m] | 0,v_mat[m,m]) + normal_lpdf(sigma2[3] |sigma2[2]* rho1, tar_sigma2) + normal_lpdf(beta[k] |0,1) + normal_lpdf(epsilon2[m,3] | 0, tau_epsilon2);
      }
    }
  }

  // 2016
  for (m in 1:M){
    for (k in 1:K){
      if (y16[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,4,k]) + normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[4] | sigma1[3]* rho1, tar_sigma1);
      } else {
        target += log(p[m,4,k]) + poisson_lpmf(to_int(y16[m,k]) | theta[m,4,k])+ normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[4] | sigma1[3]* rho1, tar_sigma1) +normal_lpdf(gamma2[m] | 0,v_mat[m,m]) + normal_lpdf(sigma2[4] |sigma2[3]* rho1, tar_sigma2) + normal_lpdf(beta[k] |0,1) + normal_lpdf(epsilon2[m,4] | 0, tau_epsilon2);
      }
    }
  }
  
  // 2017
  for (m in 1:M){
    for (k in 1:K){
      if (y17[m,k] == 0){
        target += bernoulli_lpmf(0 | p[m,5,k]) + normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[5] | sigma1[4]* rho1, tar_sigma1);
      } else {
        target += log(p[m,5,k]) + poisson_lpmf(to_int(y17[m,k]) | theta[m,5,k])+ normal_lpdf(alpha1[k] |0,1) + normal_lpdf(alpha2[k] |0,1)  + normal_lpdf(gamma1[m] | 0,v_mat[m,m]) + normal_lpdf(sigma1[5] | sigma1[4]* rho1, tar_sigma1) +normal_lpdf(gamma2[m] | 0,v_mat[m,m]) + normal_lpdf(sigma2[5] |sigma2[4]* rho1, tar_sigma2) + normal_lpdf(beta[k] |0,1) + normal_lpdf(epsilon2[m,5] | 0, tau_epsilon2);
      }
    }
  }
}

generated quantities {
  real theta_post[M,N,K]; 
  real p_post[M,N,K];

  for (m in 1:M){
    for (k in 1:K){
      for (n in 1:N){
        p_post[m,n,k] = inv_logit(col(X,k)' * alpha1 + log(pop[k,n])* col(X,k)' * alpha2 + gamma1[m] + sigma1[n]);
        theta_post[m,n,k] = exp(col(X,k)' * beta + log(pop[k,n]) + gamma2[m] + sigma2[n] + epsilon2[m,n]);
      }
    }
  }
}


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
  real tar_sigma1;
  real tar_sigma2;
  real rho1;
  real rho2;
  real<lower=0> tar_gamma1;
  real<lower=0> tar_gamma2;
  real<lower=0> tar_theta;
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
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector<lower=0> [M*N*K] p;
  vector<lower=0> [M*N*K] theta;
  vector[K] alpha1;
  vector[K] alpha2;
  vector[N*M] epilson;
  vector[K] beta;
  vector[N-1] k1;
  vector[N-1] k2;
  vector[M] tar1;
  vector[M] tar2;
}

transformed parameters {
  vector [M*N*K] log_theta = log(theta);
  vector [M*N*K] logit_p = logit(p);
  vector[N] sigma1;
  vector[N] sigma2;
  vector[M] mu;
  
  sigma1[1] = 1;
  sigma2[1] = 1;
  
  for (i in 1:M){
    mu[i] = 0;
  }
  
  for (n in 2:N) {
    sigma1[n] = rho1 * sigma1[n-1] + k1[n-1];
    sigma2[n] = rho2 * sigma2[n-1] + k2[n-1];
  }
  
  // 2013
  for (m in 1:M) {
    for (k in 1:K) {
      logit_p[(m-1)*k + k] = X[,k]' * alpha1 + log(pop[k,1]) * X[,k]' * alpha2 + tar1[m] + sigma1[1];
      log_theta[(m-1)*k + k] = log(pop[k,1]) + X[,k]' * beta + tar2[m] + sigma2[1] + epilson[m];
    }
  }

  //2014
  for (m in 1:M) {
    for (k in 1:K) {
      logit_p[(m-1)*k + k + M*K] = X[,k]' * alpha1 + log(y14[m,k]) * X[,k]' * alpha2 + tar1[m] + sigma1[2];
      log_theta[(m-1)*k + k + M*K] = log(y14[m,k]) + X[,k]' * beta + tar2[m] + sigma2[2] + epilson[M+m];
    }
  }
  
  //2015
  for (m in 1:M) {
    for (k in 1:K) {
      logit_p[(m-1)*k + k + 2*M*K] = X[,k]' * alpha1 + log(y15[m,k]) * X[,k]' * alpha2 + tar1[m] + sigma1[3];
      log_theta[(m-1)*k + k + 2*M*K] = log(y15[m,k]) + X[,k]' * beta + tar2[m] + sigma2[3] + epilson[2*M+m];
    }
  }
  
  //2016
  for (m in 1:M) {
    for (k in 1:K) {
      logit_p[(m-1)*k + k + 3*M*K] = X[,k]' * alpha1 + log(y16[m,k]) * X[,k]' * alpha2 + tar1[m] + sigma1[4];
      log_theta[(m-1)*k + k + 3*M*K] = log(y16[m,k]) + X[,k]' * beta + tar2[m] + sigma2[4] + epilson[3*M+m];
    }
  }
  
  //2017
  for (m in 1:M) {
    for (k in 1:K) {
      logit_p[(m-1)*k + k + 4*M*K] = X[,k]' * alpha1 + log(y17[m,k]) * X[,k]' * alpha2 + tar1[m] + sigma1[5];
      log_theta[(m-1)*k + k + 4*M*K] = log(y17[m,k]) + X[,k]' * beta + tar2[m] + sigma2[5] + epilson[4*M+m];
    }
  }


}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // priors
  sigma1[1] ~ normal(0, tar_sigma1*(1-rho1));
  sigma2[1] ~ normal(0, tar_sigma2*(1-rho2));
  k1 ~ normal(0, tar_sigma1);
  k2 ~ normal(0, tar_sigma2);
  tar1 ~ multi_normal(mu, tar_gamma1*C);
  tar2 ~ multi_normal(mu, tar_gamma2*C);
  epilson ~ normal(0, tar_theta);
  alpha1 ~ normal(0, 0.001);
  alpha2 ~ normal(0, 0.001);
  beta ~ normal(0, 0.001);
  // posterior
  
  //2013
  for (m in 1:M) {
    for (k in 1:K) {
      if (y13[m, k] > 0) {
        target += log(p[(m-1)*k + k]*((theta[(m-1)*k + k]^y13[m,k])*exp(-theta[(m-1)*k + k]))/(tgamma(y13[m,k]+1)*(1-exp(-theta[(m-1)*k + k]))));
      } else {
        target += log(1-p[(m-1)*k + k]);
      }
    }
  }
  
  //2014
  for (m in 1:M) {
    for (k in 1:K) {
      if (y14[m, k] > 0) {
        target += log(p[(m-1)*k + k + M*K]*((theta[(m-1)*k + k + M*K]^y14[m,k])*exp(-theta[(m-1)*k + k + M*K]))/(tgamma(y14[m,k]+1)*(1-exp(-theta[(m-1)*k + k + M*K]))));
      } else {
        target += log(1-p[(m-1)*k + k + M*K]);
      }
    }
  }
  
  //2015
  for (m in 1:M) {
    for (k in 1:K) {
      if (y15[m, k] > 0) {
        target += log(p[(m-1)*k + k + 2*M*K]*((theta[(m-1)*k + k + 2*M*K]^y15[m,k])*exp(-theta[(m-1)*k + k + 2*M*K]))/(tgamma(y15[m,k]+1)*(1-exp(-theta[(m-1)*k + k + 2*M*K]))));
      } else {
        target += log(1-p[(m-1)*k + k + 2*M*K]);
      }
    }
  }
  
  //2016
  for (m in 1:M) {
    for (k in 1:K) {
      if (y16[m, k] > 0) {
        target += log(p[(m-1)*k + k + 3*M*K]*((theta[(m-1)*k + k + 3*M*K]^y16[m,k])*exp(-theta[(m-1)*k + k + 3*M*K]))/(tgamma(y16[m,k]+1)*(1-exp(-theta[(m-1)*k + k + 3*M*K]))));
      } else {
        target += log(1-p[(m-1)*k + k + 3*M*K]);
      }
    }
  }
  
  //2017
  for (m in 1:M) {
    for (k in 1:K) {
      if (y17[m, k] > 0) {
        target += log(p[(m-1)*k + k + 4*M*K]*((theta[(m-1)*k + k + 4*M*K]^y17[m,k])*exp(-theta[(m-1)*k + k + 4*M*K]))/(tgamma(y17[m,k]+1)*(1-exp(-theta[(m-1)*k + k + 4*M*K]))));
      } else {
        target += log(1-p[(m-1)*k + k + 4*M*K]);
      }
    }
  }
}

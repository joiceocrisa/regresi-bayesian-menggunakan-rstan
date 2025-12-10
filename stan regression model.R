# ---
# Stan Regression Model
# ---

library(rstan)

stan_pois_code <-"
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> y[N];
  matrix[N,K] X;
}
parameters {
  real alpha;       // intercept
  vector[K] beta;   // slope
}
model {
  alpha ~ normal(0, 1);
  beta  ~ normal(0, 1);

  y ~ poisson_log(alpha + X * beta);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = poisson_log_lpmf(y[n] | alpha + X[n] * beta);
}
"

stan_nb_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=0> y[N];
  matrix[N,K] X;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> phi;
}
model {
  // Prior weakly-informative (standar)
  alpha ~ normal(0, 10);
  beta  ~ normal(0, 10);

  // phi > 0, prior cocok gamma/exponential
  phi   ~ gamma(2, 0.1);

  // Negative Binomial likelihood (log link)
  y ~ neg_binomial_2_log(alpha + X * beta, phi);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] =
      neg_binomial_2_log_lpmf(y[n] | alpha + X[n] * beta, phi);
}
"


stan_mix_pois_code <- "
data {
  int<lower=1> N;
  int<lower=0> y[N];
  int<lower=1> K;
  matrix[N, K] X;
}
parameters {
  simplex[2] theta;        // bobot campuran (theta[1] + theta[2] = 1)
  ordered[2] alpha;        // intercept untuk komponen 1 & 2 (diurutkan agar identifiable)
  vector[K] beta1;         // koefisien regresi komponen 1
  vector[K] beta2;         // koefisien regresi komponen 2
}
transformed parameters {
  vector[N] eta1;
  vector[N] eta2;
  for (n in 1:N) {
    eta1[n] = alpha[1] + X[n] * beta1;  // log(lambda_1)
    eta2[n] = alpha[2] + X[n] * beta2;  // log(lambda_2)
  }
}
model {
  // Priors
  alpha ~ normal(0, 10);
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  
  // Likelihood mixture
  for (n in 1:N) {
    target += log_mix(theta[1],
                      poisson_log_lpmf(y[n] | eta1[n]),
                      poisson_log_lpmf(y[n] | eta2[n]));
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = log_mix(theta[1],
                         poisson_log_lpmf(y[n] | eta1[n]),
                         poisson_log_lpmf(y[n] | eta2[n]));
  }
}
"


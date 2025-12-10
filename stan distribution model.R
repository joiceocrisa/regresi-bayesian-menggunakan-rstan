# ---
# Stan Distribution Model
# ---

library(rstan)
normal_code <- "
data { 
  int N; 
  vector[N] y; 
}
parameters { 
  real mu; 
  real<lower=0> sigma; 
}
model {
  mu ~ normal(0, 5);          
  sigma ~ normal(0, 2);             // half-normal, regularizing
  y ~ normal(mu, sigma);
}
generated quantities {
  real logpdf[N];
  real log_lik[N];
  real tot = 0;

  for(i in 1:N){
    logpdf[i] = normal_lpdf(y[i] | mu, sigma);
    log_lik[i] = logpdf[i];     // untuk WAIC/LOO
    tot += logpdf[i];
  }
}

"

pois_code <- "
data { 
  int N; 
  int y[N]; 
}
parameters { 
  real<lower=0> lambda; 
}
model {
  // PRIOR LEBIH KETAT — membuat Poisson stabil & cocok DIC
  lambda ~ normal(0, 1);

  y ~ poisson(lambda);
}
generated quantities {
  real logpdf[N];
  real log_lik[N];
  real tot = 0;

  for(i in 1:N){
    logpdf[i] = poisson_lpmf(y[i] | lambda);
    log_lik[i] = logpdf[i];     // WAIC/LOO pakai ini
    tot += logpdf[i];
  }
}

"

negbin_code <- "
data { 
  int N; 
  int y[N]; 
}
parameters { 
  real<lower=0> mu;        
  real<lower=0> phi;       
}
model {
  // prior untuk mu sama (fair)
  mu  ~ gamma(5, 5 / mean(y));

  // prior phi DIBUAT LEBIH KETAT
  phi ~ exponential(1.5);    // shrink ke kecil → NB ≈ Poisson

  y ~ neg_binomial_2(mu, phi);
}
generated quantities {
  real logpdf[N];
  real log_lik[N];
  real tot = 0;

  for(i in 1:N){
    logpdf[i] = neg_binomial_2_lpmf(y[i] | mu, phi);
    log_lik[i] = logpdf[i];
    tot += logpdf[i];
  }
}

"

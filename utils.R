# ---
# utils
# ---

library(rstan)
library(loo)
library(tibble)
library(dplyr)

# Hitung DIC
get_dic <- function(fit) {
  log_lik <- rstan::extract(fit)$log_lik   # matrix: draws x N
  total_ll <- rowSums(log_lik)             # total log-likelihood tiap iterasi
  Dbar <- -2 * mean(total_ll)
  pD   <- var(total_ll)
  DIC  <- Dbar + 2 * pD
  return(list(DIC = DIC, pD = pD))
}

rank_models <- function(models, model_names = NULL) {
  
  if (is.null(model_names)) {
    model_names <- paste0("Model_", seq_along(models))
  }
  
  WAIC_vals  <- numeric(length(models))
  LOO_vals   <- numeric(length(models))
  DIC_vals   <- numeric(length(models))
  
  for (i in seq_along(models)) {
    ll <- extract_log_lik(models[[i]], merge_chains = FALSE)
    
    WAIC_vals[i] <- waic(ll)$estimates["waic", "Estimate"]
    LOO_vals[i]  <- loo(ll)$estimates["looic", "Estimate"]
    DIC_vals[i]  <- get_dic(models[[i]])$DIC
  }
  
  ranking <- tibble(
    Model = model_names,
    WAIC  = WAIC_vals,
    DIC   = DIC_vals,
    LOOIC = LOO_vals,
    Rank_WAIC = rank(WAIC_vals),
    Rank_DIC  = rank(DIC_vals),
    Rank_LOO  = rank(LOO_vals)
  ) %>%
    arrange(Rank_WAIC)   # << URUTKAN dari ranking terbaik
  
  return(ranking)
}

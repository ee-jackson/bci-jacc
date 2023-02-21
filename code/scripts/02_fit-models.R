#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: fit-models
## Desc: fit models and save output
## Date created: 2023-02-06


# packages ----------------------------------------------------------------

library("tidyverse")
library("here")
library("rstan")
library("brms")
library("bayestestR")

# data --------------------------------------------------------------------

readRDS(here::here("data", "clean", "connect_pod_data.rds")) %>% 
  mutate(connectivity_sc = scale(connectivity),
         individ_fecundity_sc = scale(individ_fecundity),
         dbh_mm_sc = scale(dbh_mm),
         n_mature_sc = scale(n_mature),
         n_total_sc = scale(n_total)) %>% 
  mutate(n_predated = round(n_predated), n_total = round(n_total)
  ) -> model_data

# Bayesian model ----------------------------------------------------------

bprior <- prior(normal(0, 1), class = b)

full_model <-
  brm(data = model_data,
      family = binomial(link = logit),
      n_predated | trials(n_total) ~
        connectivity_sc + individ_fecundity_sc +
        connectivity_sc:individ_fecundity_sc,
      prior = bprior,
      iter = 12500,
      warmup = 500,
      chains = 4,
      cores = 4,
      seed = 9,
      file = (here::here("output", "models", "model_fit.rds")))

bayestestR::describe_prior(full_model) %>% 
  as.data.frame() %>%
  write_csv(here::here("output", "results", "describe_prior.csv"))

bayestestR::describe_posterior(full_model,
                               ci = 0.95,
                               ci_method = "HDI",
                               test = c("p_direction", "rope"),
                               centrality = "median") %>%
  as.data.frame() %>%
  write_csv(here::here("output", "results", "describe_posterior.csv"))


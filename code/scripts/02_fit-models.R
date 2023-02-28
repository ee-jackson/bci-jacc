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

options(brms.file_refit = "on_change")


# data --------------------------------------------------------------------

readRDS(here::here("data", "clean", "connect_pod_data.rds")) %>% 
  drop_na() %>%
  mutate(across(contains("n_"), ~as.integer(.)) %>% 
  mutate(connectivity_sc = scale(connectivity),
         individ_fecundity_sc = scale(individ_fecundity),
         dbh_mm_sc = scale(dbh_mm),
         n_mature_sc = scale(n_mature),
         n_total_sc = scale(n_total)) 
  ) -> model_data


# Bayesian model ----------------------------------------------------------

bprior <- prior(normal(0, 1), class = b)

model_predation <-
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
      file = (here::here("output", "models", "predation_model_fit.rds"))
      )

bayestestR::describe_prior(model_predation) %>% 
  as.data.frame() %>%
  write_csv(here::here("output", "results", "predation_model_describe_prior.csv"))

bayestestR::describe_posterior(model_predation,
                               ci = 0.95,
                               ci_method = "HDI",
                               test = c("p_direction", "rope"),
                               centrality = "median") %>%
  as.data.frame() %>%
  write_csv(here::here("output", "results", "predation_model_describe_posterior.csv"))


# Fruit set model ---------------------------------------------------------
model_fruit_set <-
  brm(data = model_data,
      family = gaussian,
      n_total_sc ~
        connectivity_sc + dbh_mm_sc,
      prior = bprior,
      iter = 12500,
      warmup = 500,
      chains = 4,
      cores = 4,
      seed = 9,
      file = (here::here("output", "models", "fruit_set_model_fit.rds")))

bayestestR::describe_posterior(model_fruit_set,
                               ci = 0.95,
                               ci_method = "HDI",
                               centrality = "median") %>% 
  as.data.frame() %>% 
  write_csv(here::here("output", "results", "fruit_set_model_describe_posterior.csv"))


# Realised fecundity model ------------------------------------------------

model_realised_fecundity <-
  brm(data = model_data,
      family = gaussian,
      n_mature_sc ~
        connectivity_sc + dbh_mm_sc,
      prior = bprior,
      iter = 12500,
      warmup = 500,
      chains = 4,
      cores = 4,
      seed = 9,
      file = (here::here("output", "models", "realised_fecundity_model_fit.rds")))

bayestestR::describe_posterior(model_realised_fecundity,
                               ci = 0.95,
                               ci_method = "HDI",
                               centrality = "median") %>% 
  as.data.frame() %>% 
  write_csv(here::here("output", "results", "realised_fecundity_model_describe_posterior.csv"))

# Fruit length model ------------------------------------------------------

read.csv(here::here("data", "raw", "jacaranda_pods.csv"),
         header = TRUE, na.strings = c("", "NA", "missing")) %>%
  filter(tree != "JACC_130") %>%
  mutate(fragment = case_when(str_detect(comments, "fragment") ~ TRUE,
                              TRUE ~ FALSE)) %>%  
  mutate(pod_half_whole = recode(pod_half_whole, "half" = 0.5,
                                 "whole" = 1)) %>%
  mutate(pod_size_mm = coalesce(pod_size_string_mm, pod_size_mm)) %>%
  mutate(pod_size_mm = as.numeric(pod_size_mm)) %>% 
  filter(fragment == FALSE) %>% 
  mutate(pod_mature = case_when(
    pod_size_mm >= 40 & str_detect(morph, "^symmetrical_locules") ~ TRUE,
    is.na(pod_size_mm) ~ NA,
    is.na(morph) ~ NA,
    TRUE ~ FALSE
  )) %>%
  filter(pod_mature == TRUE) -> mature_pods
  
mature_pods %>% 
  filter(pod_half_whole == 1) %>% 
  bind_rows(mature_pods) -> mature_halves

mature_halves %>% 
  left_join(model_data, by = c("tree" = "tree_id"), multiple = "all") %>% 
  mutate(pod_size_sc = scale(pod_size_mm)) -> length_data

model_length <-
  brm(data = length_data,
      family = gaussian,
      pod_size_sc ~
        connectivity_sc + individ_fecundity_sc + (1|tree),
      prior = bprior,
      iter = 12500,
      warmup = 500,
      chains = 4,
      cores = 4,
      seed = 9,
      file = (here::here("output", "models", "length_model_fit.rds")))

bayestestR::describe_posterior(model_length,
                               ci = 0.95,
                               ci_method = "HDI",
                               centrality = "median") %>% 
  as.data.frame() %>% 
  write_csv(here::here("output", "results", "length_model_describe_posterior.csv"))

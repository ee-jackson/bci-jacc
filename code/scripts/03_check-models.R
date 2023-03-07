#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: check-models
## Desc: perform posterior predictive checks
## Date created: 2023-03-06

# packages ----------------------------------------------------------------

library("tidyverse")
library("here")
library("brms")
library("patchwork")
library("stringr")
library("bayesplot")


# Posterior predictive checks ---------------------------------------------

file_names <- as.list(dir(path = here::here("output", "models"),
                          pattern = "*.rds", full.names = TRUE))

model_list <- lapply(file_names, readRDS)

names(model_list) <- lapply(file_names, basename)

plot_pp_check <- function(model) {
  pp_check(model, ndraws = 50) +  
    ggtitle(str_wrap(formula(model), 60)) +
    theme_classic(base_size = 5)
}

plot_list <- lapply(model_list, plot_pp_check)

wrap_plots(plot_list)

ggsave(here::here("output","figures","pp_checks.png"),
       width = 1476, height = 1476, units = "px")


# MCMC diagnostics --------------------------------------------------------

plot_mcmc_check <- function(model) {
  mcmc_trace(model, regex_pars = "b_",
             n_warmup = 500) +  
    ggtitle(str_wrap(formula(model), 70)) +
    theme_classic(base_size = 5)
}

mcmc_plot_list <- lapply(model_list, plot_mcmc_check)

wrap_plots(mcmc_plot_list, ncol = 1)

ggsave(here::here("output","figures","mcmc_checks.png"),
       width = 1476, height = 2000, units = "px")


# Check influence of prior information ------------------------------------

prior <- distribution_normal(n = 50, mean = 0, sd = 1)
plot(density(x))

ggplot() +
  geom_density(data = posteriors, aes(x = b_individ_fecundity_sc), fill = "orange") +
  geom_density(aes(x = prior))


intercept_prior <- distribution_student_t(n = 50, df = 3, ncp = 0)

ggplot() +
  geom_density(data = posteriors, aes(x = b_Intercept), fill = "orange") +
  geom_density(aes(x = intercept_prior))

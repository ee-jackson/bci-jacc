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
library("bayestestR")
library("gt")

file_names <- as.list(dir(path = here::here("output", "models"),
                          pattern = "*.rds", full.names = TRUE))

model_list <- lapply(file_names, readRDS)

names(model_list) <- lapply(file_names, basename)

# Posterior predictive checks ---------------------------------------------

plot_pp_check <- function(model, xlab) {
  pp_check(model, ndraws = 500) +
    labs(x = xlab, y = "Density") +
    theme_classic(base_size = 20)
}


# MCMC diagnostics --------------------------------------------------------

plot_mcmc_check <- function(model) {
  mcmc_trace(model, regex_pars = "b_",
             iter1 = 500,
             facet_args = list(ncol = 1)) +
    theme_classic(base_size = 15)
}

# Get posterior param estimates -------------------------------------------

get_table <- function(model) {
  bayestestR::describe_posterior(model,
                                 ci = 0.95,
                                 ci_method = "HDI",
                                 centrality = "median",
                                 test = FALSE) %>%
    mutate(across(!Rhat & !Parameter, round, 2)) %>%
    gt()
}


# fruit set model ---------------------------------------------------------

fs_t <- get_table(model_list$fruit_set_model_fit.rds)
gtsave(fs_t, here::here("output", "figures", "fruit_set_model_out.png"))
fs_t_png <- png::readPNG(here::here("output", "figures", "fruit_set_model_out.png"),
                         native = TRUE)

fs_pp <- plot_pp_check(model_list$fruit_set_model_fit.rds, 
                       xlab = "Individual fecundity")
fs_mcmc <- plot_mcmc_check(model_list$fruit_set_model_fit.rds)

((fs_pp / fs_t_png) | fs_mcmc) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))

png(
  here::here("output", "figures", "fruit_set_si.png"),
  width = 1476,
  height = 1000,
  units = "px"
)


# fruit size model --------------------------------------------------------

l_t <- get_table(model_list$length_model_fit.rds)
gtsave(l_t, here::here("output", "figures", "fruit_length_model_out.png"))
l_t_png <- png::readPNG(here::here("output", "figures", "fruit_length_model_out.png"),
                         native = TRUE)

l_pp <- plot_pp_check(model_list$length_model_fit.rds, 
                      xlab = "Fruit size (mm)")
l_mcmc <- plot_mcmc_check(model_list$length_model_fit.rds)

((l_pp / l_t_png) | l_mcmc) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))

png(
  here::here("output", "figures", "fruit_length_si.png"),
  width = 1476,
  height = 1000,
  units = "px"
)


# predation model ---------------------------------------------------------

pred_t <- get_table(model_list$predation_model_fit.rds)
gtsave(pred_t, here::here("output", "figures", "fruit_pred_model_out.png"))
pred_t_png <- png::readPNG(here::here("output", "figures", "fruit_pred_model_out.png"),
                        native = TRUE)

pred_pp <- plot_pp_check(model_list$predation_model_fit.rds, 
                         xlab = "Predated fruit count")
pred_mcmc <- plot_mcmc_check(model_list$predation_model_fit.rds)

((pred_pp / pred_t_png) | pred_mcmc) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))

png(
  here::here("output", "figures", "fruit_pred_si.png"),
  width = 1476,
  height = 1000,
  units = "px"
)


# mature fruit model ------------------------------------------------------

mat_t <- get_table(model_list$realised_fecundity_model_fit.rds)
gtsave(mat_t, here::here("output", "figures", "fruit_mature_model_out.png"))
mat_t_png <- png::readPNG(here::here("output", "figures", "fruit_mature_model_out.png"),
                           native = TRUE)

mat_pp <- plot_pp_check(model_list$realised_fecundity_model_fit.rds, 
                        xlab = "Realised fecundity")
mat_mcmc <- plot_mcmc_check(model_list$realised_fecundity_model_fit.rds)

((mat_pp / mat_t_png) | mat_mcmc) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))

png(
  here::here("output", "figures", "fruit_mature_si.png"),
  width = 1476,
  height = 1000,
  units = "px"
)

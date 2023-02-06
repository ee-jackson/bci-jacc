#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: fit-models
## Desc:
## Date created: 2023-02-06


# packages ----------------------------------------------------------------

library("here")
library("tidyverse")
library("lme4")
#library("DHARMa")
#library("broom.mixed")


# data --------------------------------------------------------------------

readRDS(here::here("data", "clean", "connect_pod_data.rds")) -> full_data


# J & C model
full_data %>%
  drop_na(n_predated) %>%
  mutate(n_predated = round(n_predated), n_total = round(n_total)) %>%
  lme4::glmer(
            formula = cbind(n_predated, (n_total - n_predated)) ~ 
              log(individ_fecundity) + log(connectivity) + 
              log(individ_fecundity):log(connectivity),
            family = binomial(link = "logit")
) -> m1

summary(m1)

#
full_data %>%
  drop_na(n_predated) %>%
  mutate(n_predated = round(n_predated), n_total = round(n_total)) %>%
  glm(
    formula = cbind(n_predated, (n_total - n_predated)) ~ 
      log(connectivity),
    family = binomial(link = "logit")
  ) -> m4

summary(m4)

#
full_data %>%
  drop_na(n_predated) %>%
  mutate(n_predated = round(n_predated), n_total = round(n_total)) %>%
  glm(
    formula = cbind(n_predated, (n_total - n_predated)) ~ 
      log(individ_fecundity),
    family = binomial(link = "logit")
  ) -> m5

summary(m5)

# effect of connectivity on individual fruit production
glm(full_data,
            formula = log(n_total) ~ 
              log(connectivity) + log(dbh_mm),
            family = gaussian(link = "identity")
) -> m2

summary(m2)

# effect of connectivity on realised fruit production
full_data %>%
  drop_na(n_mature) %>%
  glm(full_data,
            formula = log(n_mature) ~ 
              log(connectivity) + log(dbh_mm),
            family = gaussian(link = "identity")
) -> m3

summary(m3)

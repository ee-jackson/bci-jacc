#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: summarise-data
## Desc: make summary tables and figures
## Date created: 2023-02-06


# packages ----------------------------------------------------------------

library("tidyverse")
library("here")
library("modelsummary")
library("kableExtra")

# data --------------------------------------------------------------------

readRDS(here::here("data", "clean", "connect_pod_data.rds")) -> full_data
  
full_data %>%
  rowwise() %>%
  drop_na() %>%
  mutate(`Proportion of fruit predated per tree` = n_predated / n_total) %>%
  rename(
    `Estimated capsule production per tree` = individ_fecundity,
    `Mature capsules per tree` = n_mature
  ) -> new_data

tmp <- new_data[, c("Estimated capsule production per tree",
                    "Proportion of fruit predated per tree",
                    "Mature capsules per tree")]

tmp_list <- lapply(tmp, na.omit)

emptycol = function(x) " "

datasummary(
  `Estimated capsule production per tree` +
  `Proportion of fruit predated per tree` +
  `Mature capsules per tree` ~ Mean + SD + 
    Heading("Boxplot") * emptycol + 
    Heading("Histogram") * emptycol, 
  data = tmp) %>%
  kable_classic() %>%
  column_spec(column = 1, width = "10em") %>%
  column_spec(column = 4, image = spec_boxplot(tmp_list,
    width = 700, col = "#6baed6",  medlwd = 2, medcol = "red",
    height = 300, add_label = TRUE, same_lim = FALSE
  )) %>%
  column_spec(column = 5, image = spec_hist(tmp_list,
    width = 700, col = "#6baed6", border = "#6baed6",
    height = 200, same_lim = FALSE
  )) %>%
  kable_styling(font_size = 16, html_font = "arial") %>% 
  kableExtra::as_image(file = here::here("output", "figures", "summary-stats.pdf"), 
                       width = 4.92) 

# for kableExtra::as_image units are in inches
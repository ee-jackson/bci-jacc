#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: make-figures
## Desc: make all figures to be included in main text
## Date created: 2023-02-06


# Packages ----------------------------------------------------------------

library("tidyverse")
library("here")
library("tidybayes")
library("ggdist")
library("modelsummary")
library("kableExtra")

full_data <- readRDS(here::here("data", "clean", "connect_pod_data.rds"))

model_fit <- readRDS(here::here("output", "models", "model_fit.rds"))


# Summarise data ----------------------------------------------------------

full_data %>%
  rowwise() %>%
  drop_na() %>%
  mutate(`Proportion of fruit predated per tree` = n_predated / n_total) %>%
  rename(
    `Estimated capsule production per tree` = individ_fecundity,
    `Mature capsules per tree` = n_mature
  ) %>%
  select(`Estimated capsule production per tree`,
         `Proportion of fruit predated per tree`,
         `Mature capsules per tree`) -> plotting_data

emptycol <- function(x) " "

datasummary(
  `Estimated capsule production per tree` +
    `Proportion of fruit predated per tree` +
    `Mature capsules per tree` ~ Mean + SD + 
    Heading("Boxplot") * emptycol + 
    Heading("Histogram") * emptycol, 
  data = plotting_data) %>%
  kable_classic() %>%
  column_spec(column = 1, width = "10em") %>%
  column_spec(column = 4, 
              image = spec_boxplot(tmp_list,
                                   width = 700, col = "#6baed6",  medlwd = 2, 
                                   medcol = "red", height = 300, add_label = TRUE, 
                                   same_lim = FALSE)) %>%
  column_spec(column = 5, 
              image = spec_hist(tmp_list,
                                width = 700, col = "#6baed6", border = "#6baed6",
                                height = 200, same_lim = FALSE)) %>%
  kable_styling(font_size = 16, html_font = "arial") %>% 
  kableExtra::as_image(file = here::here("output", "figures", "summary_stats.pdf"), 
                       width = 4.92) 

# for kableExtra::as_image units are in inches


# Visualise model results -------------------------------------------------

png(here::here("output","figures","parameter_estimates.png"), 
    width = 1476, height = 1100, units = "px")
tidybayes::tidy_draws(model_fit) %>% 
  rename(`Individual fecundity` = b_individ_fecundity_sc, 
         Connectivity = b_connectivity_sc,
         `Individual fecundity :\nConnectivity` =`b_connectivity_sc:individ_fecundity_sc`) %>% 
  select(`Individual fecundity`, Connectivity, `Individual fecundity :\nConnectivity`) %>% 
  pivot_longer(cols = everything(), names_to = "parameter") %>%
  ggplot(aes(y = reorder(parameter, abs(value)), x = value)) +
  ggdist::stat_halfeye(aes(fill = after_stat(cut_cdf_qi(cdf, .width= c(.66, .95, 1)))),
                       stroke = 2, size = 20, linewidth = 15, 
                       shape = 21, slab_alpha = 1,
                       point_fill = "white") +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues", direction=-1) +
  labs(x = "Possible parameter values", y = "") +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 2, colour = "red") +
  theme_classic(base_size = 50) +
  theme(legend.position = "none",
        plot.margin = margin(t = 0.1, r = 0.2, b = 0.1, l = 0, "in"),
        axis.text.y = element_text(colour="black")) 
dev.off()
 
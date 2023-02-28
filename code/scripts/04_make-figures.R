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
library("ggmap")

full_data <- readRDS(here::here("data", "clean", "connect_pod_data.rds"))

full_data %>% 
  drop_na() %>%
  mutate(across(contains("n_"), ~as.integer(.)) %>% 
           mutate(connectivity_sc = scale(connectivity),
                  individ_fecundity_sc = scale(individ_fecundity),
                  dbh_mm_sc = scale(dbh_mm),
                  n_mature_sc = scale(n_mature),
                  n_total_sc = scale(n_total)) 
  ) -> model_data

readRDS(here::here("data", "clean", "all_tree_lon_lat.rds")) %>%
  mutate(focal = ifelse(grepl("CGL", tree_id), FALSE, TRUE)) -> all_trees

pred_mod <- readRDS(here::here("output", "models", "predation_model_fit.rds"))


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
tidybayes::tidy_draws(pred_mod) %>% 
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
        plot.margin = margin(t = 0.1, r = 0.4, b = 0.1, l = 0, "in"),
        axis.text.y = element_text(colour="black")) 
dev.off()

# connectivity
model_data %>%
  modelr::data_grid(connectivity_sc = modelr::seq_range(connectivity_sc, n = 51),
                    n_total = modelr::seq_range(n_total, n = 51),
                    individ_fecundity_sc = rep(mean(model_data$individ_fecundity_sc), 51)) %>% 
  mutate(n_total = as.integer(n_total)) %>% 
  add_epred_draws(pred_mod, ndraws = 500) %>%
  mutate(connectivity_us = connectivity_sc * 
           attr(model_data$connectivity_sc, 'scaled:scale') + attr(model_data$connectivity_sc, 'scaled:center')) %>% 
  ggplot() +
  stat_lineribbon(aes(x = connectivity_us, y = .epred/n_total), colour = "darkblue", alpha = 0.7) +
  geom_point(data = model_data, aes(x = connectivity, y = n_predated/n_total), 
             inherit.aes = FALSE, alpha = 0.7, shape = 16, size = 5) +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues") +
  theme_classic(base_size = 35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(legend.position = "none") +
  xlab("Connectivity") +
  ylab("Proportion of capsules predated") -> p1

# fecundity
model_data %>%
  modelr::data_grid(individ_fecundity_sc = modelr::seq_range(individ_fecundity_sc, n = 51),
                    n_total = modelr::seq_range(n_total, n = 51),
                    connectivity_sc = rep(mean(model_data$connectivity_sc), 51)) %>% 
  mutate(n_total = as.integer(n_total)) %>% 
  add_epred_draws(pred_mod, ndraws = 500) %>%
  mutate(individ_fecundity_us = individ_fecundity_sc * 
           attr(model_data$individ_fecundity_sc, 'scaled:scale') + attr(model_data$individ_fecundity_sc, 'scaled:center'),) %>% 
  ggplot() +
  stat_lineribbon(aes(x = individ_fecundity_us, y = .epred/n_total), colour = "darkblue", alpha = 0.7) +
  geom_point(data = model_data, aes(x = individ_fecundity, y = n_predated/n_total), 
             inherit.aes = FALSE, alpha = 0.7, shape = 16, size = 5) +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues") +
  theme_classic(base_size = 35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(legend.position = "none") +
  xlab("Individual fecundity") +
  ylab("") -> p2

png(here::here("output","figures","epred_draws.png"), 
    width = 1476, height = 700, units = "px")
p1 + p2
dev.off()

# Plot trees on density map -----------------------------------------------

bbox <- make_bbox(c(min(all_trees$lon) - 0.001, max(all_trees$lon) + 0.001), 
                  c(min(all_trees$lat) - 0.001, max(all_trees$lat) + 0.001))

bci_basemap <- ggmap::get_map(bbox,
                              force = TRUE, maptype = "toner-lite")

all_trees %>% 
  left_join(full_data) %>% 
  group_by(focal) %>% 
  group_split() -> split_trees

bci_basemap %>% 
  ggmap() +
  geom_point(data = all_trees, aes(lon, lat),
             size = 2, shape = 16, colour = "#6baed6", alpha = 0.5) +
  geom_point(data = split_trees[[2]], aes(lon, lat, size = connectivity),
             colour = "red", shape = 4) +
  theme_void(base_size = 7) +
  theme(legend.position = c(.99, .01),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill="transparent", colour="transparent")) -> p_map

ggsave(here::here("output","figures","map3.png"), plot = p_map,
       width = 1476, height = 1476, units = "px")


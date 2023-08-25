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
library("patchwork")

full_data <- readRDS(here::here("data", "clean", "connect_pod_data.rds"))

full_data %>% 
  drop_na() %>%
  mutate(across(contains("n_"), ~as.integer(.)) %>% 
           mutate(connectivity_sc = scale(connectivity),
                  individ_fecundity_sc = scale(individ_fecundity),
                  dbh_mm_sc = scale(dbh_mm)) 
  ) -> model_data

readRDS(here::here("data", "clean", "all_tree_lon_lat.rds")) %>%
  mutate(focal = ifelse(grepl("CGL", tree_id), FALSE, TRUE)) -> all_trees


# Summarise data ----------------------------------------------------------

full_data %>%
  rowwise() %>%
  drop_na() %>%
  mutate(`Proportion of fruits predated` = n_predated / n_total) %>%
  rename(
    `Estimated fruit set` = individ_fecundity,
    `Count of mature fruits` = n_mature
  ) %>%
  select(`Estimated fruit set`,
         `Proportion of fruits predated`,
         `Count of mature fruits`) -> plotting_data

emptycol <- function(x) " "

tmp_list <- as.list(plotting_data)

datasummary(
  `Estimated fruit set` +
    `Proportion of fruits predated` +
    `Count of mature fruits` ~ Mean + SD + 
    Heading("Boxplot") * emptycol + 
    Heading("Histogram") * emptycol, 
  data = plotting_data) %>%
  kable_classic() %>%
  column_spec(column = 1, width = "8em") %>%
  column_spec(column = 4, 
              image = spec_boxplot(tmp_list,
                                   width = 800, col = "#6baed6",  medlwd = 2, 
                                   medcol = "red", height = 400, add_label = TRUE, 
                                   same_lim = FALSE)) %>%
  column_spec(column = 5, 
              image = spec_hist(tmp_list,
                                width = 800, col = "#6baed6", border = "#6baed6",
                                height = 300, same_lim = FALSE)) %>%
  kable_styling(font_size = 25, html_font = "arial") %>% 
  kableExtra::as_image(file = here::here("output", "figures", "summary_stats.png"), 
                       width = 10) 

# for kableExtra::as_image units are in inches

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
             colour = "darkblue", shape = 4) +
  theme_void(base_size = 7) +
  theme(legend.position = c(.99, .01),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill="transparent", colour="transparent")) -> p_map

ggsave(here::here("output","figures","map.png"), plot = p_map,
       width = 1476, height = 1476, units = "px")

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
                  dbh_mm_sc = scale(dbh_mm),
                  n_mature_sc = scale(n_mature),
                  n_total_sc = scale(n_total)) 
  ) -> model_data

readRDS(here::here("data", "clean", "all_tree_lon_lat.rds")) %>%
  mutate(focal = ifelse(grepl("CGL", tree_id), FALSE, TRUE)) -> all_trees


# Summarise data ----------------------------------------------------------

full_data %>%
  rowwise() %>%
  drop_na() %>%
  mutate(`Proportion of fruits predated` = n_predated / n_total) %>%
  rename(
    `Individual fecundity` = individ_fecundity,
    `Realised fecundity` = realised_fecundity
  ) %>%
  select(`Individual fecundity`,
         `Proportion of fruits predated`,
         `Realised fecundity`) -> plotting_data

emptycol <- function(x) " "
Min <- function(x) sprintf("%.2f", min(x, na.rm = TRUE))
Max <- function(x) sprintf("%.2f", max(x, na.rm = TRUE))

tmp_list <- as.list(plotting_data)

# make inset plots
png(here::here("output", "plots","hist_indi_fec.png"),
    width = 300, height = 200, units = "px", bg = "transparent")
par(mar = c(2, 2, 1, 2))
hist(plotting_data$`Individual fecundity`, 
     col = "#6baed6", border = NA, main = NA, xlab = NA, ylab = NA)
dev.off()

png(here::here("output", "plots","hist_prop_pred.png"),
    width = 300, height = 200, units = "px", bg = "transparent")
par(mar = c(2, 2, 1, 2))
hist(plotting_data$`Proportion of fruits predated`, 
     col = "#6baed6", border = NA, main = NA, xlab = NA, ylab = NA)
dev.off()

png(here::here("output", "plots","hist_real_fec.png"),
    width = 300, height = 200, units = "px", bg = "transparent")
par(mar = c(2, 2, 1, 2))
hist(plotting_data$`Realised fecundity`, 
     col = "#6baed6", border = NA, main = NA, xlab = NA, ylab = NA)
dev.off()

png(here::here("output", "plots","box_indi_fec.png"),
    width = 300, height = 200, units = "px", bg = "transparent")
par(mar = c(2, 2, 1, 2))
boxplot(plotting_data$`Individual fecundity`, frame = FALSE,
        horizontal = TRUE, col = "#6baed6", main = NA, xlab = NA, ylab = NA)
dev.off()

png(here::here("output", "plots","box_prop_pred.png"),
    width = 300, height = 200, units = "px", bg = "transparent")
par(mar = c(2, 2, 1, 2))
boxplot(plotting_data$`Proportion of fruits predated`, frame = FALSE,
        horizontal = TRUE, col = "#6baed6", main = NA, xlab = NA, ylab = NA)
dev.off()

png(here::here("output", "plots","box_real_fec.png"),
    width = 300, height = 200, units = "px", bg = "transparent")
par(mar = c(2, 2, 1, 2))
boxplot(plotting_data$`Realised fecundity`, frame = FALSE,
        horizontal = TRUE, col = "#6baed6", main = NA, xlab = NA, ylab = NA)
dev.off()

# make table
datasummary(
  `Individual fecundity` +
    `Proportion of fruits predated` +
    `Realised fecundity` ~ Min + Max + Mean + SD +
    Heading("Boxplot") * emptycol + 
    Heading("Histogram") * emptycol, 
  data = plotting_data) %>%
  kable_classic() %>%
  column_spec(column = 1, width_min = "5em") %>%
  kableExtra::column_spec(column = 6, image = c(
    here::here("output", "plots","box_indi_fec.png"),
    here::here("output", "plots","box_prop_pred.png"),
    here::here("output", "plots","box_real_fec.png"),
    here::here("output", "plots","box_real_fec.png")
  )) %>% 
  kableExtra::column_spec(column = 7, image = c(
    here::here("output", "plots","hist_indi_fec.png"),
    here::here("output", "plots","hist_prop_pred.png"),
    here::here("output", "plots","hist_real_fec.png"),
    here::here("output", "plots","hist_real_fec.png")
               )) %>% 
  kable_styling(font_size = 30, html_font = "arial")%>% 
  kableExtra::as_image(file = here::here("output", "figures", "summary_stats.png"), 
                       width = 4.92) 

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

#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-connectivity
## Desc:
## Date created: 2023-01-23


# packages ----------------------------------------------------------------

library("tidyverse")
library("here")
library("geosphere")
library("sf")


# data --------------------------------------------------------------------

readRDS(here::here("data", "clean", "pod_data.rds")) %>% 
  mutate(crown_area_m2 = pi * crown_radius_m^2) %>% 
  select(tree_id, lon, lat, crown_area_m2) -> focal_jacc

st_read(here::here(
  "data", "maps",
  "jacc-map-garzonlopez2012", "JAC1COpointSept.shp"
)) %>% 
  filter(!st_is_empty(.)) %>%
  select(geometry) -> point_data

st_read(here::here(
  "data", "maps",
  "polygons-garzonlopez2008", "JAC1CO_pol2008.shp"
)) %>% 
  filter(!st_is_empty(.)) %>%
  select(geometry) -> polygon_data

# filter CGL data to trees that have polygon data and calculate crown area (m^2)
st_join(polygon_data, point_data, 
        join = st_within, largest = TRUE, left = FALSE) %>%  
  mutate(crown_area_m2 = as.numeric(st_area(geometry))) %>% 
  mutate(tree_id = paste0("CGL_", 1:n()),
         geometry = st_centroid(geometry)) %>% 
  sf::st_transform(crs = "+proj=longlat +datum=WGS84", allow_ballpark = FALSE) %>% 
  mutate(lon = sf::st_coordinates(.)[,1],
         lat = sf::st_coordinates(.)[,2]) -> other_jacc
  
# keep only adult trees
other_jacc %>%
  mutate(dbh_estimate_mm = (crown_area_m2 + 99.42) / 0.38) %>%
  filter(dbh_estimate_mm >= 200) %>% 
  select(- dbh_estimate_mm) -> other_jacc

# remove duplicate trees --------------------------------------------------

geosphere::distm(cbind(pull(focal_jacc, lon), pull(focal_jacc, lat)),
  cbind(pull(other_jacc, lon), pull(other_jacc, lat)),
  fun = distGeo
) -> dist_matrix

rownames(dist_matrix) <- focal_jacc$tree_id
colnames(dist_matrix) <- other_jacc$tree_id

# assume closest CGL tree within 20 m of focal tree is the same tree
as.data.frame.table(dist_matrix, responseName = "dist") %>%
  filter(dist <= 20) %>% 
  group_by(Var1) %>%
  slice(which.min(dist)) %>% 
  ungroup() -> duplicate_trees
  
other_jacc %>%
  filter(!tree_id %in% pull(duplicate_trees, Var2)) -> other_jacc_no_dupes

other_jacc_no_dupes %>% 
  st_drop_geometry() %>% 
  bind_rows(focal_jacc) -> all_jacc

# if focal tree was duplicated in CGL data, use polygon for crown area not radius
duplicate_trees %>%
  left_join(other_jacc, by = c("Var2" = "tree_id")) %>% 
  select(Var1, crown_area_m2) %>% 
  rename(tree_id = Var1) %>%
  right_join(y = all_jacc, by = "tree_id", 
             suffix = c("_poly", "_radius")) %>% 
  mutate(crown_area_m2 = ifelse(is.na(crown_area_m2_poly), 
                                crown_area_m2_radius, crown_area_m2_poly)) %>%
  select(- crown_area_m2_poly, -crown_area_m2_radius) -> all_jacc_areas
  
  
# calculate pairwise distances --------------------------------------------

calculate_dist <- function(data) {
  data %>%
    select(lon, lat) -> plot_matrix

  geosphere::distm(plot_matrix,
    fun = distGeo
  ) -> dists

  as.data.frame(dists) -> dists_df

  unlist(data$tree_id) -> colnames(dists_df)

  cbind(data, dists_df)
}

calculate_dist(all_jacc_areas) -> distance_df


# calculate connectivity --------------------------------------------------

calculate_connectivity <- function(data, id) {
  data %>%
    filter(tree_id != eval(parse(text = id))) %>%
    select(id, crown_area_m2) %>%
    mutate(x = exp(-1 / 85 * eval(parse(text = id))) * crown_area_m2) %>%
    summarise(
      tree_id = paste(id),
      connectivity = sum(x)
    )
}

focal_jacc %>%
  distinct(tree_id) %>%
  pull(tree_id) -> focal_id_list

lapply(focal_id_list,
  calculate_connectivity,
  data = distance_df) -> connectivity_dfs

connectivity_dfs %>%
  dplyr::bind_rows() -> all_connectivity_dfs


# join and save -----------------------------------------------------------

readRDS(here::here("data", "clean", "pod_data.rds")) %>%
  left_join(all_connectivity_dfs, by = "tree_id") %>%
  saveRDS(here::here("data", "clean", "connect_pod_data.rds"))

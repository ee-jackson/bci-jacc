#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: calculate-connectivity
## Desc: 
## Date created: 2023-01-23


# packages ----------------------------------------------------------------

library("tidyverse") 
library("here") 
library("geosphere")

# data --------------------------------------------------------------------

readRDS(here::here("data", "clean", "pod_data.rds")) %>%
  select(tree_id, lon, lat) -> focal_jacc

plotKML::readGPX(here::here("data", "maps",
                            "jacc-map-garzonlopez2012", "jac1co_map.gpx")) %>%
  map_df( ~ .) %>%
  select(lon, lat) %>%
  mutate(tree_id = paste0("CGL_", 1:n())) -> other_jacc


# remove duplicate trees --------------------------------------------------

geosphere::distm(cbind(pull(focal_jacc, lon), pull(focal_jacc, lat)), 
                 cbind(pull(other_jacc, lon), pull(other_jacc, lat)), 
                 fun = distGeo) -> dist_matrix

rownames(dist_matrix) <- focal_jacc$tree_id
colnames(dist_matrix) <- other_jacc$tree_id

# assume closest tree within 30 m from focal tree is the same tree
as.data.frame.table(dist_matrix, responseName = "dist") %>% 
  filter(dist <= 30) %>%
  group_by(Var1) %>%
  slice(which.min(dist)) %>% 
  ungroup() %>%
  pull(Var2) %>%
  unique() -> duplicate_trees

rbind(other_jacc, focal_jacc) %>% 
  filter(!tree_id %in% duplicate_trees) -> all_jacc


# calculate pairwise distances --------------------------------------------

calculate_dist <- function (data) {
  
  data %>%
    select(lat, lon) -> plot_matrix
  
  rdist::pdist(plot_matrix[,c("lat", "lon")], 
               metric = "euclidean") -> dists
  
  as.data.frame(dists) -> dists_df
  
  unlist(data$tree_id) -> colnames(dists_df) 
  
  cbind(data, dists_df)
  
}

distance_df <- calculate_dist(all_jacc)


# calculate connectivity --------------------------------------------------

focal_jacc %>%
  distinct(tree_id) %>%
  pull(tree_id) -> focal_tree_id_list

calculate_connectivity <- function (data, id) {
  data %>%
    select(id) %>%
    mutate(x = exp(- eval(parse(text = id)) ) ) %>%
    summarise(tree_id = paste(id),
              connectivity = sum(x))
}

connectivity_dfs <- lapply(focal_tree_id_list,
                           calculate_connectivity, data = distance_df)

connectivity_dfs %>%
  dplyr::bind_rows() -> all_connectivity_dfs


# join and save -----------------------------------------------------------

readRDS(here::here("data", "clean", "pod_data.rds")) %>%
  left_join(all_connectivity_dfs) %>%
  saveRDS(here::here("data", "clean", "connect_pod_data.rds"))

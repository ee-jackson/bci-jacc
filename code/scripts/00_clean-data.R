#!/usr/bin/env Rscript

## Author: E E Jackson, eleanor.elizabeth.j@gmail.com
## Script: clean-data
## Desc: get the raw data ready for analysis
## Date created: 2023-01-19


# packages ----------------------------------------------------------------

library("tidyverse") 
library("here") 

# data --------------------------------------------------------------------

read.csv(here::here("data", "raw", "jacaranda_pods.csv"),
         header = TRUE, na.strings = c("", "NA", "missing")) %>%
  filter(tree != "JACC_130") %>%
  mutate(fragment = case_when(str_detect(comments, "fragment") ~ TRUE,
                              TRUE ~ FALSE)) %>%  
  mutate(pod_half_whole = recode(pod_half_whole, "half" = 0.5,
                                 "whole" = 1)) %>%
  mutate(pod_size_mm = coalesce(pod_size_string_mm, pod_size_mm)) %>%
  mutate(pod_size_mm = as.numeric(pod_size_mm)) -> pod_data

# we found no capsules at JACC_130 but there is still a data entry with dbh

pod_data %>%
  select(tree, date, dbh_mm,crown_radius_m) %>%
  rename(tree_id = tree) %>%
  distinct() -> tree_data
# add crown area column


# count immature pods -----------------------------------------------------

pod_data %>%
  filter(fragment == FALSE) %>% 
  mutate(pod_immature = case_when(
    pod_size_mm < 40 & str_detect(morph, "^symmetrical_locules") ~ TRUE,
    is.na(pod_size_mm) ~ NA,
    is.na(morph) ~ NA,
    TRUE ~ FALSE
  )) %>%
  filter(pod_immature == TRUE) %>%
  group_by(tree) %>%
  summarise(n_immature = sum(pod_half_whole, na.rm = TRUE), .groups = "drop") -> immature_pods


# count mature pods -------------------------------------------------------

pod_data %>%
  filter(fragment == FALSE) %>% 
  mutate(pod_mature = case_when(
    pod_size_mm >= 40 & str_detect(morph, "^symmetrical_locules") ~ TRUE,
    is.na(pod_size_mm) ~ NA,
    is.na(morph) ~ NA,
    TRUE ~ FALSE
  )) %>%
  filter(pod_mature == TRUE) %>%
  group_by(tree) %>%
  summarise(n_mature = sum(pod_half_whole, na.rm = TRUE), .groups = "drop") -> mature_pods


# count predated pods -----------------------------------------------------

pod_data %>%
  filter(fragment == FALSE) %>% 
  mutate(pod_predated = case_when(
    morph == "asymmetrical_locules" |
      morph == "single_locule" |
      morph == "small_knobbly" ~ TRUE,
    is.na(morph) ~ NA,
    TRUE ~ FALSE
  )) %>% 
  filter(pod_predated == TRUE) %>%
  group_by(tree) %>%
  summarise(n_predated = sum(pod_half_whole, na.rm = TRUE), .groups = "drop") -> predated_pods


# count total pods --------------------------------------------------------

pod_data %>%
  filter(fragment == FALSE) %>% 
  group_by(tree) %>%
  summarise(n_total = sum(pod_half_whole, na.rm = TRUE), .groups = "drop") -> total_pods


# estimate fruit set ------------------------------------------------------

pod_data %>%
  filter(fragment == FALSE) %>% 
  select(tree, crown_radius_m) %>%
  distinct() %>%
  mutate(crown_area = pi * crown_radius_m^2) %>%
  mutate(n_quadrats = case_when(
    str_detect(crown_radius_m, "\\d{0,2}(\\.0{1})") ~ floor(crown_radius_m),
         TRUE ~ ceiling(crown_radius_m)
    )) %>%
  mutate(n_quadrats = n_quadrats * 4) %>%
  left_join(total_pods) %>%
  rowwise() %>%
  mutate(mean_fruit_set = n_total/n_quadrats) %>%
  mutate(individ_fecundity = round(mean_fruit_set * crown_area)) %>%
  select(tree, individ_fecundity) -> fruit_set


# combine groups ----------------------------------------------------------

full_join(immature_pods, mature_pods, by = "tree") %>%
  full_join(predated_pods, by = "tree") %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
  full_join(total_pods, by = "tree") %>%
  full_join(fruit_set, by = "tree") -> grouped_pod_data


# get tree co-ordinates ---------------------------------------------------

plotKML::readGPX(here::here("data", "raw", "jacaranda_waypoints.gpx")) %>%
  map_df( ~ .) %>%
  rename(tree_id = name) %>%
  arrange(tree_id) %>%
  select(tree_id, lon, lat) %>%
  filter(!grepl("NS", tree_id)) %>%
  mutate(tree_id = gsub("EJ", "", tree_id) ) %>%
  mutate(tree_id = gsub("^(JACC)", "\\1_\\2", tree_id)) %>%
  right_join(tree_data) -> tree_data_coords


# output ------------------------------------------------------------------

full_join(tree_data_coords, grouped_pod_data, by = c("tree_id" = "tree")) %>%
  saveRDS(here::here("data", "clean", "pod_data.rds"))

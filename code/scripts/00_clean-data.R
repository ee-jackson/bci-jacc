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

pod_data %>%
  filter(fragment == FALSE) %>% 
  group_by(tree) %>%
  summarise(n_total_pods = sum(pod_half_whole, na.rm = TRUE), .groups = "drop") -> total_pods
  
full_join(immature_pods, mature_pods, by = "tree") %>%
  full_join(predated_pods, by = "tree") %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
  full_join(total_pods, by = "tree") -> grouped_pod_data

full_join(tree_data, grouped_pod_data, by = c("tree_id" = "tree"))
# add individual fecundity
library("tidybayes")
library("tidyverse")
library("patchwork")

readRDS(here::here("data", "clean", "connect_pod_data.rds")) %>% 
  drop_na() %>%
  mutate(across(contains("n_"), ~as.integer(.)) %>% 
           mutate(connectivity_sc = scale(connectivity),
                  individ_fecundity_sc = scale(individ_fecundity),
                  dbh_mm_sc = scale(dbh_mm),
                  n_mature_sc = scale(n_mature),
                  n_total_sc = scale(n_total)) 
  ) -> model_data

read.csv(here::here("data", "raw", "jacaranda_pods.csv"),
         header = TRUE, na.strings = c("", "NA", "missing")) %>%
  filter(tree != "JACC_130") %>%
  mutate(fragment = case_when(str_detect(comments, "fragment") ~ TRUE,
                              TRUE ~ FALSE)) %>%  
  mutate(pod_half_whole = recode(pod_half_whole, "half" = 0.5,
                                 "whole" = 1)) %>%
  mutate(pod_size_mm = coalesce(pod_size_string_mm, pod_size_mm)) %>%
  mutate(pod_size_mm = as.numeric(pod_size_mm)) %>% 
  filter(fragment == FALSE) %>% 
  mutate(pod_mature = case_when(
    pod_size_mm >= 55 & str_detect(morph, "^symmetrical_locules") ~ TRUE,
    is.na(pod_size_mm) ~ NA,
    is.na(morph) ~ NA,
    TRUE ~ FALSE
  )) %>%
  filter(pod_mature == TRUE) -> mature_pods

mature_pods %>% 
  filter(pod_half_whole == 1) %>% 
  bind_rows(mature_pods) -> mature_halves

mature_halves %>% 
  left_join(model_data, by = c("tree" = "tree_id"), multiple = "all") %>% 
  mutate(pod_size_sc = scale(pod_size_mm)) -> length_data

# get models
file_names <- as.list(dir(path = here::here("output", "models"),
                          pattern = "*.rds", full.names = TRUE))

model_list <- lapply(file_names, readRDS)

names(model_list) <- lapply(file_names, basename)

model_data %>%
  modelr::data_grid(connectivity_sc = modelr::seq_range(connectivity_sc, n = 51),
                    n_total = modelr::seq_range(n_total, n = 51),
                    individ_fecundity_sc = rep(mean(model_data$individ_fecundity_sc), 51)) %>% 
  mutate(n_total = as.integer(n_total)) %>% 
  add_epred_draws(model_list$predation_model_fit.rds, ndraws = 500) %>%
  mutate(connectivity_us = connectivity_sc * 
           attr(model_data$connectivity_sc, 'scaled:scale') + attr(model_data$connectivity_sc, 'scaled:center')) %>% 
  ggplot() +
  geom_point(data = model_data, aes(x = connectivity, y = n_predated/n_total, size = n_total), 
             inherit.aes = FALSE, alpha = 0.7, shape = 16) +
  stat_lineribbon(aes(x = connectivity_us, y = .epred/n_total), 
                  colour = "darkblue", alpha = 0.4) +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues") +
  theme_classic(base_size = 35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(legend.position = "none") +
  xlab("Connectivity") +
  ylab("Proportion predated") -> p1

model_data %>%
  modelr::data_grid(individ_fecundity_sc = modelr::seq_range(individ_fecundity_sc, n = 51),
                    n_total = modelr::seq_range(n_total, n = 51),
                    connectivity_sc = rep(mean(model_data$connectivity_sc), 51)) %>% 
  mutate(n_total = as.integer(n_total)) %>% 
  add_epred_draws(model_list$predation_model_fit.rds, ndraws = 500) %>%
  mutate(individ_fecundity_us = individ_fecundity_sc * 
           attr(model_data$individ_fecundity_sc, 'scaled:scale') + attr(model_data$individ_fecundity_sc, 'scaled:center')) %>% 
  ggplot() +
  geom_point(data = model_data, aes(x = individ_fecundity, y = n_predated/n_total, size = n_total), 
             inherit.aes = FALSE, alpha = 0.7, shape = 16) +
  stat_lineribbon(aes(x = individ_fecundity_us, y = .epred/n_total), 
                  colour = "darkblue", alpha = 0.4) +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues") +
  theme_classic(base_size = 35) +
  scale_x_continuous(expand = c(0.01, 0.5)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(legend.position = "none") +
  xlab("Individual fecundity") +
  ylab("Proportion predated") -> p2

model_data %>%
  modelr::data_grid(connectivity_sc = modelr::seq_range(connectivity_sc, n = 51)) %>% 
  add_epred_draws(model_list$fruit_set_model_fit.rds, ndraws = 500) %>%
  mutate(connectivity_us = connectivity_sc * 
           attr(model_data$connectivity_sc, 'scaled:scale') + attr(model_data$connectivity_sc, 'scaled:center')) %>% 
  ggplot() +
  geom_point(data = model_data, aes(x = connectivity, y = n_total), 
             inherit.aes = FALSE, alpha = 0.6, shape = 16, size = 5) +
  stat_lineribbon(aes(x = connectivity_us, y = .epred), 
                  colour = "darkblue", alpha = 0.4) +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues") +
  theme_classic(base_size = 35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.005, 5)) +
  theme(legend.position = "none") +
  xlab("Connectivity") +
  ylab("Individual fecundity") -> p3

model_data %>%
  modelr::data_grid(connectivity_sc = modelr::seq_range(connectivity_sc, n = 51)) %>% 
  add_epred_draws(model_list$realised_fecundity_model_fit.rds, ndraws = 500) %>%
  mutate(connectivity_us = connectivity_sc * 
           attr(model_data$connectivity_sc, 'scaled:scale') + attr(model_data$connectivity_sc, 'scaled:center')) %>% 
  ggplot() +
  geom_point(data = model_data, aes(x = connectivity, y = n_mature), 
             inherit.aes = FALSE, alpha = 0.6, shape = 16, size = 5) +
  stat_lineribbon(aes(x = connectivity_us, y = .epred), 
                  colour = "darkblue", alpha = 0.4) +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues") +
  theme_classic(base_size = 35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.005, 5)) +
  theme(legend.position = "none") +
  xlab("Connectivity") +
  ylab("Realised fecundity") -> p5


length_data %>%
  modelr::data_grid(connectivity_sc = modelr::seq_range(connectivity_sc, n = 51),
                    individ_fecundity_sc = rep(mean(model_data$individ_fecundity_sc), 51)) %>% 
  add_epred_draws(model_list$length_model_fit.rds, ndraws = 500, re_formula = NA) %>%
  mutate(individ_fecundity_us = individ_fecundity_sc * 
           attr(model_data$individ_fecundity_sc, 'scaled:scale') + attr(model_data$individ_fecundity_sc, 'scaled:center'),
         connectivity_us = connectivity_sc * 
           attr(model_data$connectivity_sc, 'scaled:scale') + attr(model_data$connectivity_sc, 'scaled:center')) %>% 
  ggplot() +
  geom_jitter(data = length_data, aes(x = connectivity, y = pod_size_mm), 
             inherit.aes = FALSE, alpha = 0.25, shape = 16, size = 2, width = 1.5) +
  stat_summary(data = length_data, aes(x = connectivity, y = pod_size_mm, group = tree), fun = median, geom = "point",
               size = 5, stroke = 1, shape = 21,
               fill = "black", colour = "darkblue" ) +
  stat_lineribbon(aes(x = connectivity_us, y = .epred), 
                  colour = "darkblue", alpha = 0.4) +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues") +
  theme_classic(base_size = 35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(legend.position = "none") +
  xlab("Connectivity") +
  ylab("Fruit size (mm)") -> p7

length_data %>%
  modelr::data_grid(individ_fecundity_sc = modelr::seq_range(individ_fecundity_sc, n = 51),
                    connectivity_sc = rep(mean(model_data$connectivity_sc), 51)) %>% 
  add_epred_draws(model_list$length_model_fit.rds, ndraws = 500, re_formula = NA) %>%
  mutate(individ_fecundity_us = individ_fecundity_sc * 
           attr(model_data$individ_fecundity_sc, 'scaled:scale') + attr(model_data$individ_fecundity_sc, 'scaled:center'),
         connectivity_us = connectivity_sc * 
           attr(model_data$connectivity_sc, 'scaled:scale') + attr(model_data$connectivity_sc, 'scaled:center')) %>% 
  ggplot() +
  geom_jitter(data = length_data, aes(x = individ_fecundity, y = pod_size_mm), 
              inherit.aes = FALSE, alpha = 0.25, shape = 16, size = 2, width = 10) +
  stat_summary(data = length_data, aes(x = individ_fecundity, y = pod_size_mm, group = tree), fun = median, geom = "point",
               size = 5, stroke = 1, shape = 21,
               fill = "black", colour = "darkblue" ) +
  stat_lineribbon(aes(x = individ_fecundity_us, y = .epred), 
                  colour = "darkblue", alpha = 0.4) +
  scale_fill_brewer(na.translate = FALSE, palette = "Blues") +
  theme_classic(base_size = 35) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  theme(legend.position = "none") +
  xlab("Individual fecundity") +
  ylab("Fruit size (mm)") -> p8

png(
  here::here("output", "figures", "all_models_new.png"),
  width = 1476,
  height = 1800,
  units = "px",
  type = "cairo"
)
( (p2 + p1 ) / (p3 + p5) / (p8 + p7) ) +
  plot_annotation(tag_levels = "a")
dev.off()

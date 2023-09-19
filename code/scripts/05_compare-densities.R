
# Packages ----------------------------------------------------------------

library("tidyverse")
library("here")
library("sf")
library("ggtext")

readRDS(here::here("data", "clean", "all_tree_lon_lat.rds")) %>%
  mutate(focal = ifelse(grepl("CGL", tree_id), FALSE, TRUE)) -> all_jacc

read_sf(here::here("data", "maps", "BCI_Plot_50ha", 
                   "BCI_Plot_50ha.shp")) %>%
  st_transform(plot_50ha, crs = st_crs(4326)) -> plot_50ha

# calculate pairwise distances --------------------------------------------

calculate_dist <- function(data) {
  data %>%
    select(lon, lat) -> plot_matrix
  
  geosphere::distm(plot_matrix,
                   fun = geosphere::distGeo
  ) -> dists
  
  as.data.frame(dists) -> dists_df
  
  unlist(data$tree_id) -> colnames(dists_df)
  
  cbind(data, dists_df)
}

calculate_dist(all_jacc) -> distance_df


# calculate connectivity --------------------------------------------------

calculate_connectivity <- function(data, id) {
  
  data[data$tree_id != eval(id),] %>%
    select(id, crown_area_m2) %>%
    mutate(x = exp(-1 / 120 * eval(parse(text = id))) * crown_area_m2^(0.5)) %>%
    summarise(
      tree_id = paste(id),
      connectivity = sum(x)
    )
}

all_jacc %>%
  distinct(tree_id) %>%
  pull(tree_id) -> focal_id_list

lapply(focal_id_list,
       calculate_connectivity,
       data = distance_df) -> connectivity_dfs

connectivity_dfs %>%
  dplyr::bind_rows() -> all_connectivity_dfs

read_sf(here::here("data", "maps",
                   "jacc-map-garzonlopez2012", "JAC1COpointSept.shp")) %>%
  st_transform(crs = st_crs(4326)) %>%
  sf::st_intersection(plot_50ha) %>% 
  mutate(lon = sf::st_coordinates(.)[,1],
         lat = sf::st_coordinates(.)[,2],
         fdp = TRUE) %>% 
  st_drop_geometry() %>% 
  select(lat, lon, fdp) -> jacc_fdp

all_jacc %>% 
  left_join(jacc_fdp, by = c("lat", "lon")) %>% 
  left_join(all_connectivity_dfs, by = "tree_id") %>% 
  mutate(group = case_when(
    focal == TRUE ~ "<i> J. copaia</i> sampled in this study",
    fdp == TRUE ~ "<i> J. copaia</i> in the forest dynamics plot",
    .default = "Other <i> J. copaia</i> from across BCI"
  )) -> plotting_data

plotting_data %>%
  ggplot(aes(x = connectivity, fill = group)) +
  geom_histogram(position = "identity",
                 alpha = 0.4,
                 bins = 50) +
  geom_step(aes(y=..count.., colour=group),
            stat = "bin", bins = 50, direction = "mid",
            show.legend = FALSE) +
  scale_color_manual(values = c("#0072B2", "#E69F00","#009E73")) +
  scale_fill_manual(values = c("#0072B2",  "#E69F00", "#009E73")) +
  theme_classic(base_size = 10) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Connectivity to conspecifics") +
  ylab("n trees") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 60, 10)) +
  theme_classic(base_size = 10) +
  theme(legend.title = element_blank(), 
        legend.text = element_markdown(),
        legend.position = c(.75, .80))

ggsave(here::here("output","figures","jacc_hist.png"),
       width = 1476, height = 1000, units = "px")

plotting_data %>%
  ggplot(aes(x = connectivity, fill = group)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("#0072B2",  "#E69F00", "#009E73"))+
  theme_classic(base_size = 7) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Connectivity to conspecifics") +
  ylab("Density") +
  theme(legend.title = element_blank(), 
        legend.text = element_markdown(),
        legend.position = "bottom")

ggsave(here::here("output","figures","jacc_dens.png"),
       width = 1476, height = 1000, units = "px")

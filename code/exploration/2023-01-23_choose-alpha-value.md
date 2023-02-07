Choose an alpha value
================
Eleanor Jackson
23 January, 2023

To determine the buffer radius to use for our connectivity calculation
we will fit models with varying values of *r* and compare their AIC
values. *r* is the average migration distance of our seed predator and
*r* = 1 / $\alpha$

**Connectivity**  

$$C_{i} = \sum exp(-\ \alpha \ dist_{ji})$$

``` r
library("tidyverse")
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5      ✓ purrr   0.3.4 
    ## ✓ tibble  3.1.6      ✓ dplyr   1.0.10
    ## ✓ tidyr   1.2.0      ✓ stringr 1.4.0 
    ## ✓ readr   2.0.2      ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library("geosphere")
```

## Get data

``` r
readRDS(here::here("data", "clean", "pod_data.rds")) %>%
  select(tree_id, lon, lat) -> focal_jacc

plotKML::readGPX(here::here("data", "maps",
                            "jacc-map-garzonlopez2012", "jac1co_map.gpx")) %>%
  map_df( ~ .) %>%
  select(lon, lat) %>%
  mutate(tree_id = paste0("CGL_", 1:n())) -> other_jacc
```

    ## Registered S3 methods overwritten by 'stars':
    ##   method             from
    ##   st_bbox.SpatRaster sf  
    ##   st_crs.SpatRaster  sf

## Remove duplicated trees

Trees mapped in our ‘on the ground’ data could be the same individuals
as those in Carol’s aerial maps.

``` r
geosphere::distm(cbind(pull(focal_jacc, lon), pull(focal_jacc, lat)), 
                 cbind(pull(other_jacc, lon), pull(other_jacc, lat)), 
                 fun = distGeo) -> dist_matrix

rownames(dist_matrix) <- focal_jacc$tree_id
colnames(dist_matrix) <- other_jacc$tree_id

# assume closest tree within 30 m of focal tree is the same tree
as.data.frame.table(dist_matrix, responseName = "dist") %>% 
  filter(dist <= 30) %>%
  group_by(Var1) %>%
  slice(which.min(dist)) %>% 
  ungroup() %>%
  pull(Var2) %>%
  unique() -> duplicate_trees

rbind(other_jacc, focal_jacc) %>% 
  filter(!tree_id %in% duplicate_trees) -> all_jacc
```

## Calculate pairwise distances between trees

``` r
calculate_dist <- function (data) {
  
  data %>%
    select(lon, lat) -> plot_matrix
  
  geosphere::distm(plot_matrix, 
               fun = distGeo) -> dists
  
  as.data.frame(dists) -> dists_df
  
  unlist(data$tree_id) -> colnames(dists_df) 
  
  cbind(data, dists_df)
  
}

distance_df <- calculate_dist(all_jacc)
```

## Calculate connectivity using a range of radii

``` r
focal_jacc %>%
  distinct(tree_id) %>%
  pull(tree_id) -> focal_tree_id_list

tibble(radii = seq(25, 125, 5)) -> radii_list

lapply(radii_list, function(x) 1 / x ) -> alpha_list

merge(focal_tree_id_list, alpha_list) %>%
  rename(id = x, alpha = radii) -> alpha_tree_id

calculate_connectivity <- function (data, id, alpha) {
  data %>%
    filter(tree_id != eval(parse(text = id))) %>%
    select(id) %>%
    mutate(x = exp(- alpha * eval(parse(text = id)) ) ) %>%
    summarise(connectivity = sum(x),
              tree_id = paste(id),
              alpha = eval(parse(text = paste0(alpha))))
}

connectivity_dfs <- pmap(alpha_tree_id, .f = calculate_connectivity, 
                         data = distance_df)
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(id)` instead of `id` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

``` r
connectivity_dfs %>%
  dplyr::bind_rows() -> all_connectivity_dfs

readRDS(here::here("data", "clean", "pod_data.rds")) %>% 
  drop_na() %>% 
  mutate_at(c("n_total", "n_predated"), round) %>%
  left_join(all_connectivity_dfs, by = "tree_id") %>% 
  group_by(alpha) %>%
  nest() %>%
  ungroup() %>%
  pull(data) -> connectivity_dfs_alpha
```

## Fit models and look at AIC

``` r
fit_glmms <- function (data, predictor) {
  lme4::glmer(
    data,
    formula = cbind(n_predated, (n_total - n_predated)) ~ 
      eval(parse(text = paste0("log(", predictor, ")", sep = ""))) 
    + (1 | dbh_mm),
    family = "binomial"
  )
  
}

models <- lapply(connectivity_dfs_alpha, fit_glmms, predictor = "connectivity")

lapply(models, AIC) -> AIC_list

tibble(radius = unlist(radii_list), AIC = unlist(AIC_list)) -> AIC_df

AIC_df %>%
  ggplot(aes(x = radius, y = AIC)) + geom_point()
```

![](figures/2023-01-23_choose-alpha-value/unnamed-chunk-4-1.png)<!-- -->

A radius of **75** has the lowest AIC.

I just want to check that we’re not missing something by not looking at
bigger radii. We don’t know how far our seed predator or pollinators
travel and a big part of this study was to look for J-C effects at
larger scales.

``` r
tibble(radii = seq(100, 15000, 100)) -> radii_list

lapply(radii_list, function(x) 1 / x ) -> alpha_list

merge(focal_tree_id_list, alpha_list) %>%
  rename(id = x, alpha = radii) -> alpha_tree_id

connectivity_dfs <- pmap(alpha_tree_id, .f = calculate_connectivity, 
                         data = distance_df)

connectivity_dfs %>%
  dplyr::bind_rows() -> all_connectivity_dfs

readRDS(here::here("data", "clean", "pod_data.rds")) %>% 
  drop_na() %>% 
  mutate_at(c("n_total", "n_predated"), round) %>%
  left_join(all_connectivity_dfs, by = "tree_id") %>% 
  group_by(alpha) %>%
  nest() %>%
  ungroup() %>%
  pull(data) -> connectivity_dfs_alpha

models <- lapply(connectivity_dfs_alpha, fit_glmms, predictor = "connectivity")

lapply(models, AIC) -> AIC_list

tibble(radius = unlist(radii_list), AIC = unlist(AIC_list)) %>%
  ggplot(aes(x = radius, y = AIC)) + geom_point()
```

![](figures/2023-01-23_choose-alpha-value/unnamed-chunk-5-1.png)<!-- -->

OK nope! Lots of model convergence warnings too.

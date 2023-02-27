Explore capsule sizes
================
eleanorjackson
24 February, 2023

When collecting data in the field we didn’t classify capsules as
immature or mature, but instead recorded their size and whether they
were open or not.

From this information we should be able to classify pods as immature or
mature. Immature pods do not dehisce in the canopy but might have opened
on the forest floor.

``` r
library("tidyverse")
library("here")
```

``` r
read.csv(here::here("data", "raw", "jacaranda_pods.csv"),
         header = TRUE, na.strings = c("", "NA", "missing")) %>%
  filter(tree != "JACC_130") %>%
  mutate(fragment = case_when(str_detect(comments, "fragment") ~ TRUE,
                              TRUE ~ FALSE)) %>%  
  filter(fragment == FALSE) %>% 
  mutate(pod_size_mm = coalesce(pod_size_string_mm, pod_size_mm)) %>%
  mutate(pod_size_mm = as.numeric(pod_size_mm)) -> pod_data
```

``` r
pod_data %>%
  drop_na(pod_size_mm) %>% 
  ggplot(aes(x = pod_size_mm, fill = pod_open_closed)) +
  geom_density(alpha = 0.6) +
  xlab("Capsule length (mm)") +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = 40, linetype = 2)
```

![](figures/2023-02-24_explore-capsule-lengths/open-closed-size-1.png)<!-- -->

That’s a really nice bimodal distribution. 40 mm looks like a good cut
off for immature fruits. There’s a small number of pods that are open
but are comparatively small (the first blue hump). These are probably
immature capsules that opened on the forest floor.

Let’s take a quick look at capsule sizes between the different
categories: mature, immature and predated.

``` r
pod_data %>%
  mutate(pod_category = case_when(
    pod_size_mm < 40 & str_detect(morph, "^symmetrical_locules") ~ "immature",
    pod_size_mm >= 40 & str_detect(morph, "^symmetrical_locules") ~ "mature",
    morph == "asymmetrical_locules" |
      morph == "single_locule" |
      morph == "small_knobbly" ~ "predated",
    TRUE ~ NA
  )) %>% 
  filter(!is.na(pod_category)) -> pod_groups

pod_groups %>% 
  ggplot(aes(x = pod_size_mm, fill = pod_category)) +
  geom_density(alpha = 0.6) +
  xlab("Capsule length (mm)") +
  theme(legend.title = element_blank())
```

![](figures/2023-02-24_explore-capsule-lengths/grouped-pod-size-1.png)<!-- -->

The predated capsules also have a bimodal distribution - this’ll be the
two different morphotypes that we found: “small-knobbly” ones that
looked like they’d been predated when they were still young, and
“asymmetric” ones which were bigger but had one shrunken locule.

Let’s pull some numbers out for the manuscript.

``` r
pod_groups %>% 
  group_by(pod_category) %>% 
  summarise(median = median(pod_size_mm), mean = mean(pod_size_mm))
```

    ## # A tibble: 3 × 3
    ##   pod_category median  mean
    ##   <chr>         <dbl> <dbl>
    ## 1 immature       25.9  25.8
    ## 2 mature         79.7  80.1
    ## 3 predated       28.1  35.2

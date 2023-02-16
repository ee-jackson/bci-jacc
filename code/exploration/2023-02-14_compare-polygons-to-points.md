Compare polygons to points
================
Eleanor Jackson
14 February, 2023

Carol has given us polygons for the Jacaranda so that we can calculate
crown area. But she has said that the polygon data is older than the
point data - lets see how many trees have been recruited or have died
in-between these two time points.

``` r
library("tidyverse")
library("here")
library("sf")
```

``` r
st_read(here::here(
  "data", "maps",
  "jacc-map-garzonlopez2012", "JAC1COpointSept.shp"
)) %>% 
  filter(!st_is_empty(.)) -> point_data
```

    ## Reading layer `JAC1COpointSept' from data source 
    ##   `/Users/eleanorjackson/Library/CloudStorage/OneDrive-UniversityofReading/bci-jacc/data/maps/jacc-map-garzonlopez2012/JAC1COpointSept.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 977 features and 26 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: 624357.3 ymin: 1009741 xmax: 629582.1 ymax: 1014470
    ## Projected CRS: NAD27 / UTM zone 17N

``` r
st_read(here::here(
  "data", "maps",
  "polygons-garzonlopez2008", "JAC1CO_pol2008.shp"
)) %>% 
  filter(!st_is_empty(.)) -> polygon_data
```

    ## Reading layer `JAC1CO_pol2008' from data source 
    ##   `/Users/eleanorjackson/Library/CloudStorage/OneDrive-UniversityofReading/bci-jacc/data/maps/polygons-garzonlopez2008/JAC1CO_pol2008.shp' 
    ##   using driver `ESRI Shapefile'
    ## replacing null geometries with empty geometries
    ## Simple feature collection with 990 features and 5 fields (with 4 geometries empty)
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 624350.9 ymin: 1009737 xmax: 629584.4 ymax: 1014483
    ## Projected CRS: NAD27 / UTM zone 17N

The polygon data has 986 rows and the point data has 977 rows - 9 fewer.

``` r
st_join(polygon_data, point_data, join = st_within, largest = TRUE, left = FALSE) %>% 
  nrow()
```

    ## [1] 960

There are 960 trees that are in both datasets. To calculate the number
of trees that have died we can do 986 - 960 = 26. The number that have
been recruited is 977 - 960 = 17.

We don’t want to include any dead trees in our dataset so will have to
use the crown sizes of the trees that appear in both datasets. For the
trees in the point data for which we don’t have polygons, we could use
an average crown size?

Extract, transform, and load data
================
Beau Larkin
2022-12-08

- <a href="#description" id="toc-description">Description</a>
  - <a href="#pre-cleaning" id="toc-pre-cleaning">Pre-cleaning</a>
- <a href="#package-and-library-installation"
  id="toc-package-and-library-installation">Package and library
  installation</a>
- <a href="#data" id="toc-data">Data</a>
  - <a href="#data-and-variable-views" id="toc-data-and-variable-views">Data
    and variable views</a>
    - <a href="#terpene_meta"
      id="toc-terpene_meta"><code>terpene_meta</code></a>
    - <a href="#terpene" id="toc-terpene"><code>terpene</code></a>
    - <a href="#tree_height" id="toc-tree_height"><code>tree_height</code></a>
    - <a href="#tree_meta" id="toc-tree_meta"><code>tree_meta</code></a>
    - <a href="#tree_rust_response"
      id="toc-tree_rust_response"><code>tree_rust_response</code></a>

# Description

This script pulls data from separate files in the **Database** folder,
adds a few new variables, and checks data for consistency. This script
will be called from other analysis scripts so that the ETL process isn’t
repeated, making the analysis scripts a little shorter.

## Pre-cleaning

``` r
rm(list = ls()) 
```

# Package and library installation

Note that messages and code are often hidden in this notebook for
brevity.

``` r
# Package and library installation
packages_needed <- c("tidyverse", "magrittr", "knitr")
packages_installed <- packages_needed %in% rownames(installed.packages())
```

``` r
if (any(! packages_installed))
  install.packages(packages_needed[! packages_installed])
```

``` r
for (i in 1:length(packages_needed)) {
  library(packages_needed[i], character.only = T)
}
```

# Data

- Data are extracted from the associated directory and loaded into a
  list.
- The `resistance_class` variable is synthesized and added to several
  tables

``` r
{
  fpaths <-
    list.files(paste0(getwd(), "/Database"), full.names = TRUE)
  fnames <-
    gsub(".csv", "", list.files(paste0(getwd(), "/Database"), full.names = FALSE))
  data <-
    lapply(fpaths, function(x)
      read_csv(x, show_col_types = FALSE))
  names(data) <- fnames
  # Assign resistance classes to families in `terpene`, `terpene_meta` and `tree_meta`:
  assign_resistance <- function(x) {
    x %<>% 
      mutate(resistance_class = case_when(
        family %in% c("ENDO-155", "ENDO-157", "ENDO-158") ~ "QDR",
        family %in% c("ENDO-159", "ENDO-160") ~ "susceptible",
        TRUE ~ "MGR"))
  }
  data$terpene %<>% assign_resistance()
  data$terpene_meta %<>% assign_resistance()
  data$tree_meta %<>% assign_resistance()
}
```

## Data and variable views

### `terpene_meta`

Experimental design and metadata on trees that had terpenes extracted.
Useful for an **env** explanatory file to go along with a distance
matrix produced from `terpene`.

``` r
data$terpene_meta %>% glimpse()
```

    ## Rows: 720
    ## Columns: 7
    ## $ tree_ID          <dbl> 1002, 1003, 1007, 1009, 1012, 1016, 1019, 1021, 1023,…
    ## $ year             <dbl> 2019, 2019, 2019, 2019, 2019, 2019, 2019, 2019, 2019,…
    ## $ treatment        <chr> "FFE", "FFE", "FFE", "FFE", "FFE", "FFE", "FFE", "FFE…
    ## $ assessment       <chr> "pre_rust", "pre_rust", "pre_rust", "pre_rust", "pre_…
    ## $ block            <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,…
    ## $ family           <chr> "ENDO-159", "ENDO-159", "ENDO-156", "ENDO-156", "ENDO…
    ## $ resistance_class <chr> "susceptible", "susceptible", "MGR", "MGR", "QDR", "Q…

### `terpene`

Granular terpene masses extracted from a subset of experimental
seedlings. Also includes experimental design and metadata for
convenience. Normally, only dry weights of terpenes are used in
analysis, but wet weights are also included here.

``` r
data$terpene %>% glimpse()
```

    ## Rows: 31,680
    ## Columns: 11
    ## $ tree_ID          <dbl> 1002, 1002, 1002, 1002, 1002, 1002, 1003, 1003, 1003,…
    ## $ year             <dbl> 2019, 2019, 2019, 2019, 2019, 2019, 2019, 2019, 2019,…
    ## $ treatment        <chr> "FFE", "FFE", "FFE", "FFE", "FFE", "FFE", "FFE", "FFE…
    ## $ assessment       <chr> "pre_rust", "pre_rust", "pre_rust", "pre_rust", "pre_…
    ## $ block            <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ family           <chr> "ENDO-159", "ENDO-159", "ENDO-159", "ENDO-159", "ENDO…
    ## $ class            <chr> "diterpene", "diterpene", "diterpene", "diterpene", "…
    ## $ compound         <chr> "dehydroabietic", "levopiramic", "sandaracopiramic", …
    ## $ mass_type        <chr> "dw", "dw", "dw", "dw", "dw", "dw", "dw", "dw", "dw",…
    ## $ mass             <dbl> 0.421, 8.626, 4.115, 0.324, 2.005, 1.233, 0.244, 10.8…
    ## $ resistance_class <chr> "susceptible", "susceptible", "susceptible", "suscept…

### `tree_height`

Only applicable height measurements are included here to simplify use of
height data in correlations.

``` r
data$tree_height %>% glimpse()
```

    ## Rows: 1,023
    ## Columns: 4
    ## $ tree_ID <dbl> 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 10…
    ## $ ht6     <chr> "25", "43", "33.5", "31", "44", "32.5", "38.5", "40.5", "62", …
    ## $ ht5     <dbl> 20.0, 29.0, 30.0, 23.5, 32.0, 21.0, 29.0, 24.5, 33.0, 22.5, 23…
    ## $ ht1     <dbl> 9.5, 15.3, 13.4, 10.2, 15.0, 16.3, 15.1, 15.8, 16.0, 16.5, 11.…

### `tree_meta`

Experimental design and metadata on all trees (not just those that had
terpenes extracted).

``` r
data$tree_meta %>% glimpse()
```

    ## Rows: 1,023
    ## Columns: 6
    ## $ tree_ID          <dbl> 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009,…
    ## $ family           <chr> "ENDO-159", "ENDO-159", "ENDO-159", "ENDO-159", "ENDO…
    ## $ block            <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ endo_trt         <chr> "FFE", "FFE", "FFE", "FFE", "FFE", "FFE", "FFE", "FFE…
    ## $ rust_trt         <chr> "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes…
    ## $ resistance_class <chr> "susceptible", "susceptible", "susceptible", "suscept…

### `tree_rust_response`

Disease response traits on all seedlings.

``` r
data$tree_rust_response %>% glimpse()
```

    ## Rows: 509
    ## Columns: 28
    ## $ tree_ID   <dbl> 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, …
    ## $ inoc_dens <dbl> 3400, 3400, 3400, 3400, 3400, 3400, 3400, 3400, 3400, 3400, …
    ## $ ht6       <dbl> 25.0, 43.0, 33.5, NA, NA, NA, 38.5, NA, 62.0, 39.5, 33.0, NA…
    ## $ dm6       <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, …
    ## $ sv6       <dbl> 6, 7, 7, 9, 9, 9, 7, 9, 0, 1, 3, 9, 3, 7, 9, 5, 9, 9, 3, 8, …
    ## $ alive6    <chr> "Yes", "Yes", "Yes", "No", "No", "No", "Yes", "No", "Yes", "…
    ## $ vig6      <dbl> 1, 6, 2, 3, 3, 3, 6, 3, 1, 1, 1, 3, 1, 6, 3, 1, 3, 3, 1, 6, …
    ## $ ht5       <dbl> 20.0, 29.0, 30.0, 23.5, 32.0, 21.0, 29.0, 24.5, 33.0, 22.5, …
    ## $ dm5       <dbl> 4, 4, 4, 4, 4, 4, 4, 0, 0, 4, 4, 4, 0, 4, 4, 0, 4, 4, 0, 4, …
    ## $ sv5       <dbl> 3, 5, 5, 5, 6, 1, 3, 0, 0, 1, 1, 2, 0, 3, 7, 0, 6, 6, 0, 6, …
    ## $ vig5      <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, …
    ## $ bi5       <dbl> 2, 3, 3, 3, 10, 0, 2, 0, 0, 0, 1, 4, 0, 5, 4, 0, 3, 3, 0, 6,…
    ## $ nc5       <dbl> 2, 5, 6, 3, 10, 2, 3, 0, 0, 1, 1, 4, 0, 7, 5, 0, 4, 3, 0, 7,…
    ## $ pbr5      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ br5       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ ss5       <dbl> 2, 5, 6, 3, 10, 2, 3, 0, 0, 1, 1, 4, 0, 7, 5, 0, 4, 3, 0, 7,…
    ## $ dm4       <dbl> 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 4, 0, 4, …
    ## $ sv4       <dbl> 3, 5, 3, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 1, 6, 0, 6, …
    ## $ ss4       <dbl> 2, 6, 4, 2, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 3, 3, 0, 6,…
    ## $ dm3       <dbl> 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 4, 0, 4, …
    ## $ sv3       <dbl> 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 5, 0, 3, …
    ## $ vig3      <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ bi3       <dbl> 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 5, …
    ## $ nc3       <dbl> 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 2, 0, 5, …
    ## $ pbr3      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ br3       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ ss3       <dbl> 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 2, 0, 5, …
    ## $ ht1       <dbl> 9.5, 15.3, 13.4, 10.2, 15.0, 16.3, 15.1, 15.8, 16.0, 16.5, 1…

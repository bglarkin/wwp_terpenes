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

This step clears the environment whenever data are loaded from
individual scripts, preventing collisions or use of ghost objects.

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
data[1]
```

    ## $terpene_meta
    ## # A tibble: 720 × 7
    ##    tree_ID  year treatment assessment block family   resistance_class
    ##      <dbl> <dbl> <chr>     <chr>      <dbl> <chr>    <chr>           
    ##  1    1002  2019 FFE       pre_rust       1 ENDO-159 susceptible     
    ##  2    1003  2019 FFE       pre_rust       1 ENDO-159 susceptible     
    ##  3    1007  2019 FFE       pre_rust       1 ENDO-156 MGR             
    ##  4    1009  2019 FFE       pre_rust       1 ENDO-156 MGR             
    ##  5    1012  2019 FFE       pre_rust       1 ENDO-158 QDR             
    ##  6    1016  2019 FFE       pre_rust       1 ENDO-158 QDR             
    ##  7    1019  2019 FFE       pre_rust       1 ENDO-157 QDR             
    ##  8    1021  2019 FFE       pre_rust       1 ENDO-157 QDR             
    ##  9    1023  2019 FFE       pre_rust       1 ENDO-155 QDR             
    ## 10    1028  2019 FFE       pre_rust       1 ENDO-155 QDR             
    ## # … with 710 more rows

### `terpene`

Granular terpene masses extracted from a subset of experimental
seedlings. Also includes experimental design and metadata for
convenience. Normally, only dry weights of terpenes are used in
analysis, but wet weights are also included here.

``` r
data[2]
```

    ## $terpene
    ## # A tibble: 31,680 × 11
    ##    tree_ID  year treatment assessment block family  class compo…¹ mass_…²   mass
    ##      <dbl> <dbl> <chr>     <chr>      <dbl> <chr>   <chr> <chr>   <chr>    <dbl>
    ##  1    1002  2019 FFE       pre_rust       1 ENDO-1… dite… dehydr… dw       0.421
    ##  2    1002  2019 FFE       pre_rust       1 ENDO-1… dite… levopi… dw       8.63 
    ##  3    1002  2019 FFE       pre_rust       1 ENDO-1… dite… sandar… dw       4.12 
    ##  4    1002  2019 FFE       pre_rust       1 ENDO-1… dite… neoabi… dw       0.324
    ##  5    1002  2019 FFE       pre_rust       1 ENDO-1… dite… palust… dw       2.00 
    ##  6    1002  2019 FFE       pre_rust       1 ENDO-1… dite… abietic dw       1.23 
    ##  7    1003  2019 FFE       pre_rust       1 ENDO-1… dite… dehydr… dw       0.244
    ##  8    1003  2019 FFE       pre_rust       1 ENDO-1… dite… levopi… dw      10.8  
    ##  9    1003  2019 FFE       pre_rust       1 ENDO-1… dite… sandar… dw       2.54 
    ## 10    1003  2019 FFE       pre_rust       1 ENDO-1… dite… neoabi… dw       0.191
    ## # … with 31,670 more rows, 1 more variable: resistance_class <chr>, and
    ## #   abbreviated variable names ¹​compound, ²​mass_type

### `tree_height`

Only applicable height measurements are included here to simplify use of
height data in correlations.

``` r
data[3]
```

    ## $tree_height
    ## # A tibble: 1,023 × 4
    ##    tree_ID ht6     ht5   ht1
    ##      <dbl> <chr> <dbl> <dbl>
    ##  1    1001 25     20     9.5
    ##  2    1002 43     29    15.3
    ##  3    1003 33.5   30    13.4
    ##  4    1004 31     23.5  10.2
    ##  5    1005 44     32    15  
    ##  6    1006 32.5   21    16.3
    ##  7    1007 38.5   29    15.1
    ##  8    1008 40.5   24.5  15.8
    ##  9    1009 62     33    16  
    ## 10    1010 39.5   22.5  16.5
    ## # … with 1,013 more rows

### `tree_meta`

Experimental design and metadata on all trees (not just those that had
terpenes extracted).

``` r
data[4]
```

    ## $tree_meta
    ## # A tibble: 1,023 × 6
    ##    tree_ID family   block endo_trt rust_trt resistance_class
    ##      <dbl> <chr>    <dbl> <chr>    <chr>    <chr>           
    ##  1    1001 ENDO-159     1 FFE      Yes      susceptible     
    ##  2    1002 ENDO-159     1 FFE      Yes      susceptible     
    ##  3    1003 ENDO-159     1 FFE      Yes      susceptible     
    ##  4    1004 ENDO-159     1 FFE      Yes      susceptible     
    ##  5    1005 ENDO-159     1 FFE      Yes      susceptible     
    ##  6    1006 ENDO-156     1 FFE      Yes      MGR             
    ##  7    1007 ENDO-156     1 FFE      Yes      MGR             
    ##  8    1008 ENDO-156     1 FFE      Yes      MGR             
    ##  9    1009 ENDO-156     1 FFE      Yes      MGR             
    ## 10    1010 ENDO-156     1 FFE      Yes      MGR             
    ## # … with 1,013 more rows

### `tree_rust_response`

Disease response traits on all seedlings.

``` r
data[5]
```

    ## $tree_rust_response
    ## # A tibble: 509 × 28
    ##    tree_ID inoc_d…¹   ht6   dm6   sv6 alive6  vig6   ht5   dm5   sv5  vig5   bi5
    ##      <dbl>    <dbl> <dbl> <dbl> <dbl> <chr>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1    1001     3400  25       4     6 Yes        1  20       4     3     1     2
    ##  2    1002     3400  43       4     7 Yes        6  29       4     5     1     3
    ##  3    1003     3400  33.5     4     7 Yes        2  30       4     5     1     3
    ##  4    1004     3400  NA       4     9 No         3  23.5     4     5     1     3
    ##  5    1005     3400  NA       4     9 No         3  32       4     6     1    10
    ##  6    1006     3400  NA       4     9 No         3  21       4     1     1     0
    ##  7    1007     3400  38.5     4     7 Yes        6  29       4     3     1     2
    ##  8    1008     3400  NA       4     9 No         3  24.5     0     0     1     0
    ##  9    1009     3400  62       0     0 Yes        1  33       0     0     1     0
    ## 10    1010     3400  39.5     4     1 Yes        1  22.5     4     1     1     0
    ## # … with 499 more rows, 16 more variables: nc5 <dbl>, pbr5 <dbl>, br5 <dbl>,
    ## #   ss5 <dbl>, dm4 <dbl>, sv4 <dbl>, ss4 <dbl>, dm3 <dbl>, sv3 <dbl>,
    ## #   vig3 <dbl>, bi3 <dbl>, nc3 <dbl>, pbr3 <dbl>, br3 <dbl>, ss3 <dbl>,
    ## #   ht1 <dbl>, and abbreviated variable name ¹​inoc_dens

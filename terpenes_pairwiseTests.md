Pairwise tests of terpene composition
================
Beau Larkin

Last updated: 16 December, 2022

- <a href="#description" id="toc-description">Description</a>
- <a href="#package-and-library-installation"
  id="toc-package-and-library-installation">Package and library
  installation</a>
- <a href="#data" id="toc-data">Data</a>
- <a href="#functions" id="toc-functions">Functions</a>
  - <a href="#pairwise-comparisons-treatments-vs-controls-in-all-groups"
    id="toc-pairwise-comparisons-treatments-vs-controls-in-all-groups">Pairwise
    comparisons: treatments vs. controls in all groups</a>
    - <a href="#terms" id="toc-terms">Terms</a>
  - <a
    href="#pairwise-comparisons-induced-vs-control-in-resistance-classes-and-treatments"
    id="toc-pairwise-comparisons-induced-vs-control-in-resistance-classes-and-treatments">Pairwise
    comparisons: induced vs. control in resistance classes and
    treatments</a>
    - <a href="#terms-1" id="toc-terms-1">Terms</a>
- <a href="#results" id="toc-results">Results</a>
  - <a href="#pairwise-treatment-comparisons"
    id="toc-pairwise-treatment-comparisons">Pairwise treatment
    comparisons</a>
  - <a href="#pairwise-induced-vs-controls-comparisons"
    id="toc-pairwise-induced-vs-controls-comparisons">Pairwise induced
    vs. controls comparisons</a>

# Description

Permutation tests are used to test pairwise multivariate differences in
terpene composition within resistance classes and assessments, and
between treatments. The tests are permuted within experimental
greenhouse blocks.
[Benjamini-Hochberg](doi:10.1111/j.2517-6161.1995.tb02031.x) corrections
are used to adjust p values for multiple tests.

# Package and library installation

Note that messages and code are often hidden in this notebook for
brevity.

``` r
# Package and library installation
packages_needed <- c("tidyverse", "knitr", "vegan")
packages_installed <-
  packages_needed %in% rownames(installed.packages())
```

``` r
if (any(!packages_installed))
  install.packages(packages_needed[!packages_installed])
```

``` r
for (i in 1:length(packages_needed)) {
  library(packages_needed[i], character.only = T)
}
```

# Data

See
[data_etl.md](https://github.com/bglarkin/wwp_terpenes/blob/main/data_etl.md)
for more description of the source data. Header views of each data table
are presented here. Names only provided here.

``` r
source("data_etl.R")
```

``` r
sapply(data, function(x)
  head(x, 2))
```

    ## $terpene_meta
    ## # A tibble: 2 × 7
    ##   tree_ID  year treatment assessment block family   resistance_class
    ##     <dbl> <dbl> <chr>     <chr>      <dbl> <chr>    <chr>           
    ## 1    1002  2019 FFE       pre_rust       1 ENDO-159 susceptible     
    ## 2    1003  2019 FFE       pre_rust       1 ENDO-159 susceptible     
    ## 
    ## $terpene
    ## # A tibble: 2 × 11
    ##   tree_ID  year treat…¹ asses…² block family class compo…³ mass_…⁴  mass resis…⁵
    ##     <dbl> <dbl> <chr>   <chr>   <dbl> <chr>  <chr> <chr>   <chr>   <dbl> <chr>  
    ## 1    1002  2019 FFE     pre_ru…     1 ENDO-… dite… dehydr… dw      0.421 suscep…
    ## 2    1002  2019 FFE     pre_ru…     1 ENDO-… dite… levopi… dw      8.63  suscep…
    ## # … with abbreviated variable names ¹​treatment, ²​assessment, ³​compound,
    ## #   ⁴​mass_type, ⁵​resistance_class
    ## 
    ## $tree_height
    ## # A tibble: 2 × 4
    ##   tree_ID ht6     ht5   ht1
    ##     <dbl> <chr> <dbl> <dbl>
    ## 1    1001 25       20   9.5
    ## 2    1002 43       29  15.3
    ## 
    ## $tree_meta
    ## # A tibble: 2 × 6
    ##   tree_ID family   block endo_trt rust_trt resistance_class
    ##     <dbl> <chr>    <dbl> <chr>    <chr>    <chr>           
    ## 1    1001 ENDO-159     1 FFE      Yes      susceptible     
    ## 2    1002 ENDO-159     1 FFE      Yes      susceptible     
    ## 
    ## $tree_rust_response
    ## # A tibble: 2 × 28
    ##   tree_ID inoc_dens   ht6   dm6   sv6 alive6  vig6   ht5   dm5   sv5  vig5   bi5
    ##     <dbl>     <dbl> <dbl> <dbl> <dbl> <chr>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1    1001      3400    25     4     6 Yes        1    20     4     3     1     2
    ## 2    1002      3400    43     4     7 Yes        6    29     4     5     1     3
    ## # … with 16 more variables: nc5 <dbl>, pbr5 <dbl>, br5 <dbl>, ss5 <dbl>,
    ## #   dm4 <dbl>, sv4 <dbl>, ss4 <dbl>, dm3 <dbl>, sv3 <dbl>, vig3 <dbl>,
    ## #   bi3 <dbl>, nc3 <dbl>, pbr3 <dbl>, br3 <dbl>, ss3 <dbl>, ht1 <dbl>
    ## # ℹ Use `colnames()` to see all variable names

# Functions

## Pairwise comparisons: treatments vs. controls in all groups

This function produces group comparisons between fungal controls and
treatments for each pairwise comparison among assessments and resistance
classes. A permutation test is used in each comparison, using a blocks
as strata in the model, and a correction is made to p-values due to the
number of comparisons being done.

### Terms

- permutations = 1999
- data standardization: standardize columns (scale x to zero mean and
  unit variance)
- distance metric: euclidean

``` r
pairwise_trt <- function(c, t, p = 1999) {
  out <- NULL
  df <- data$terpene %>%
    mutate(assess_class_grp = paste(assessment, resistance_class, sep = "-")) %>%
    filter(mass_type == "dw",
           assess_class_grp == c,
           treatment %in% c("Control", t)) %>%
    select(tree_ID, treatment, assessment, block, compound, mass) %>%
    pivot_wider(
      names_from = compound,
      values_from = mass,
      values_fill = 0
    )
  X <- data.frame(df %>% select(-treatment,-assessment,-block),
                  row.names = 1)
  Y <- data.frame(df %>% select(tree_ID, treatment, block),
                  row.names = 1)
  set.seed(146)
  result <-
    data.frame(
      assess_class_grp = c,
      comparison = paste(sort(unique(df$treatment)), collapse = "_"),
      adonis2(
        scale(X) ~ treatment,
        data = Y,
        permutations = p,
        method = "euclidean",
        strata = Y$block
      )
    )[1,]
  out <- rbind(out, result)
  return(out)
}
```

## Pairwise comparisons: induced vs. control in resistance classes and treatments

This function provides group comparisons between induced and control
seedlings (**rust_inoc** vs. **rust_ctrl**) for each pairwise comparison
among resistance classes and treatments. A permutation test is used for
each comparison, and no strata are used due to the incomplete design.
P-values are corrected due to the number of comparisons being made.

### Terms

- permutations = 1999
- data standardization: standardize columns (scale x to zero mean and
  unit variance)
- distance metric: euclidean

``` r
pairwise_rust <- function(rc, t, p=1999) {
  out <- NULL
  df <- data$terpene %>%
    filter(mass_type == "dw",
           assessment != "pre_rust",
           resistance_class == rc,
           treatment == t) %>%
    select(tree_ID, treatment, assessment, compound, mass) %>%
    pivot_wider(
      names_from = compound,
      values_from = mass,
      values_fill = 0
    )
  X <- data.frame(df %>% select(-treatment,-assessment),
                  row.names = 1)
  Y <- data.frame(df %>% select(tree_ID, treatment, assessment),
                  row.names = 1)
  set.seed(146)
  result <-
    data.frame(
      resistance_class = rc,
      treatment = t,
      comparison = paste(sort(unique(df$assessment)), collapse = "-"),
      adonis2(
        scale(X) ~ assessment,
        data = Y,
        permutations = p,
        method = "euclidean"
      )
    )[1,]
  out <- rbind(out, result)
  return(out)
}
```

# Results

## Pairwise treatment comparisons

A dummy variable combining assessments and resistance_classes is needed
to reduce the number of nested calls to `lapply()`. The variable is
`assess_class_grp`:

``` r
assess_class_grp <-
  data$terpene %>%
  filter(mass_type == "dw") %>%
  mutate(class_assessment = paste(assessment, resistance_class, sep = "-")) %>%
  pull(class_assessment) %>%
  unique()
```

See output table below.

``` r
pairwise_trt_apply <-
  list(
    lapply(assess_class_grp, pairwise_trt, t = "EMF"),
    lapply(assess_class_grp, pairwise_trt, t = "FFE"),
    lapply(assess_class_grp, pairwise_trt, t = "FFE+EMF")
  )
```

``` r
pairwise_trt_tab <-
  pairwise_trt_apply %>%
  bind_rows() %>%
  mutate(
    p_val = `Pr..F.`,
    p_val_adj = round(p.adjust(p_val, "BH"), 4),
    sig_05 = case_when(p_val_adj <= 0.05 ~ "*", TRUE ~ "")
  ) %>%
  rownames_to_column(var = "delete") %>%
  select(-delete,-`Pr..F.`,-Df,-SumOfSqs) %>%
  separate(assess_class_grp, c("assessment", "resistance_class"), sep = "-") %>%
  arrange(assessment, resistance_class, comparison)
```

| assessment | resistance_class | comparison      |        R2 |         F |  p_val | p_val_adj | sig_05 |
|:-----------|:-----------------|:----------------|----------:|----------:|-------:|----------:|:-------|
| pre_rust   | MGR              | Control_EMF     | 0.0607611 |  1.423220 | 0.1250 |    0.1298 |        |
| pre_rust   | MGR              | Control_FFE     | 0.0713347 |  1.689912 | 0.0725 |    0.0851 |        |
| pre_rust   | MGR              | Control_FFE+EMF | 0.0617460 |  1.447809 | 0.0825 |    0.0928 |        |
| pre_rust   | QDR              | Control_EMF     | 0.0241587 |  1.732975 | 0.0905 |    0.0977 |        |
| pre_rust   | QDR              | Control_FFE     | 0.0329685 |  2.386472 | 0.0250 |    0.0355 | \*     |
| pre_rust   | QDR              | Control_FFE+EMF | 0.0272552 |  1.961318 | 0.0685 |    0.0841 |        |
| pre_rust   | susceptible      | Control_EMF     | 0.0461969 |  2.227982 | 0.0360 |    0.0486 | \*     |
| pre_rust   | susceptible      | Control_FFE     | 0.0387840 |  1.856051 | 0.0605 |    0.0778 |        |
| pre_rust   | susceptible      | Control_FFE+EMF | 0.0216899 |  1.019854 | 0.4040 |    0.4040 |        |
| rust_ctrl  | MGR              | Control_EMF     | 0.2743502 |  6.805353 | 0.0005 |    0.0011 | \*     |
| rust_ctrl  | MGR              | Control_FFE     | 0.2036794 |  4.603962 | 0.0005 |    0.0011 | \*     |
| rust_ctrl  | MGR              | Control_FFE+EMF | 0.1682921 |  3.642213 | 0.0005 |    0.0011 | \*     |
| rust_ctrl  | QDR              | Control_EMF     | 0.2136643 | 15.488123 | 0.0005 |    0.0011 | \*     |
| rust_ctrl  | QDR              | Control_FFE     | 0.1908055 | 13.676218 | 0.0005 |    0.0011 | \*     |
| rust_ctrl  | QDR              | Control_FFE+EMF | 0.0665173 |  4.132917 | 0.0010 |    0.0017 | \*     |
| rust_ctrl  | susceptible      | Control_EMF     | 0.2891032 | 15.860283 | 0.0005 |    0.0011 | \*     |
| rust_ctrl  | susceptible      | Control_FFE     | 0.2051627 |  9.808526 | 0.0005 |    0.0011 | \*     |
| rust_ctrl  | susceptible      | Control_FFE+EMF | 0.0725307 |  2.971709 | 0.0050 |    0.0079 | \*     |
| rust_inoc  | MGR              | Control_EMF     | 0.1625376 |  3.493502 | 0.0010 |    0.0017 | \*     |
| rust_inoc  | MGR              | Control_FFE     | 0.1492610 |  2.982626 | 0.0010 |    0.0017 | \*     |
| rust_inoc  | MGR              | Control_FFE+EMF | 0.1321500 |  2.740913 | 0.0080 |    0.0120 | \*     |
| rust_inoc  | QDR              | Control_EMF     | 0.0922763 |  5.896096 | 0.0005 |    0.0011 | \*     |
| rust_inoc  | QDR              | Control_FFE     | 0.0885755 |  5.733830 | 0.0005 |    0.0011 | \*     |
| rust_inoc  | QDR              | Control_FFE+EMF | 0.0762155 |  4.785206 | 0.0005 |    0.0011 | \*     |
| rust_inoc  | susceptible      | Control_EMF     | 0.1265003 |  5.503165 | 0.0005 |    0.0011 | \*     |
| rust_inoc  | susceptible      | Control_FFE     | 0.0965334 |  4.060214 | 0.0005 |    0.0011 | \*     |
| rust_inoc  | susceptible      | Control_FFE+EMF | 0.0965256 |  4.059851 | 0.0010 |    0.0017 | \*     |

Pairwise comparisons of terpene composition done by permutation
(n=1999); p-values corrected by the Benjamini-Hochberg method.

## Pairwise induced vs. controls comparisons

The variable `resistance_classes` is used to pass these strings to the
function via `lapply()`:

``` r
resistance_classes <- c("QDR", "susceptible", "MGR")
```

See output table below.

``` r
pairwise_rust_apply <- 
  list(
    lapply(resistance_classes, pairwise_rust, t="EMF"),
    lapply(resistance_classes, pairwise_rust, t="FFE"),
    lapply(resistance_classes, pairwise_rust, t="FFE+EMF"),
    lapply(resistance_classes, pairwise_rust, t="Control")
  )
```

``` r
pairwise_rust_tab <-
  pairwise_rust_apply %>%
  bind_rows() %>%
  mutate(
    p_val = `Pr..F.`,
    p_val_adj = round(p.adjust(p_val, "BH"), 4),
    sig_05 = case_when(p_val_adj <= 0.05 ~ "*", TRUE ~ "")
  ) %>%
  rownames_to_column(var = "delete") %>%
  select(-delete,-`Pr..F.`,-Df,-SumOfSqs)
```

| resistance_class | treatment | comparison          |        R2 |         F |  p_val | p_val_adj | sig_05 |
|:-----------------|:----------|:--------------------|----------:|----------:|-------:|----------:|:-------|
| QDR              | EMF       | rust_ctrl-rust_inoc | 0.0477863 |  2.860511 | 0.0225 |    0.0245 | \*     |
| susceptible      | EMF       | rust_ctrl-rust_inoc | 0.0731100 |  3.076188 | 0.0115 |    0.0138 | \*     |
| MGR              | EMF       | rust_ctrl-rust_inoc | 0.1121037 |  2.272638 | 0.0315 |    0.0315 | \*     |
| QDR              | FFE       | rust_ctrl-rust_inoc | 0.1101377 |  7.302397 | 0.0005 |    0.0008 | \*     |
| susceptible      | FFE       | rust_ctrl-rust_inoc | 0.1035787 |  4.390784 | 0.0005 |    0.0008 | \*     |
| MGR              | FFE       | rust_ctrl-rust_inoc | 0.1433767 |  2.845363 | 0.0060 |    0.0080 | \*     |
| QDR              | FFE+EMF   | rust_ctrl-rust_inoc | 0.1072669 |  6.969028 | 0.0005 |    0.0008 | \*     |
| susceptible      | FFE+EMF   | rust_ctrl-rust_inoc | 0.1656554 |  7.544731 | 0.0005 |    0.0008 | \*     |
| MGR              | FFE+EMF   | rust_ctrl-rust_inoc | 0.1723180 |  3.747484 | 0.0005 |    0.0008 | \*     |
| QDR              | Control   | rust_ctrl-rust_inoc | 0.1814528 | 12.857248 | 0.0005 |    0.0008 | \*     |
| susceptible      | Control   | rust_ctrl-rust_inoc | 0.1887296 |  8.840117 | 0.0005 |    0.0008 | \*     |
| MGR              | Control   | rust_ctrl-rust_inoc | 0.2218512 |  5.131824 | 0.0005 |    0.0008 | \*     |

Pairwise comparisons of terpene composition done by permutation
(n=1999); p-values corrected by the Benjamini-Hochberg method.

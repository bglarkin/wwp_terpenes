Pairwise tests of terpene composition
================
Beau Larkin
2022-12-09

- <a href="#description" id="toc-description">Description</a>
- <a href="#package-and-library-installation"
  id="toc-package-and-library-installation">Package and library
  installation</a>
- <a href="#data" id="toc-data">Data</a>
- <a href="#pairwise-comparisions-function"
  id="toc-pairwise-comparisions-function">Pairwise comparisions
  function</a>
- <a href="#results" id="toc-results">Results</a>

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
packages_needed <- c("tidyverse", "knitr", "vegan", "colorspace")
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

``` r
# Load ggplot styles and themes from text file
source("gg_style.txt")
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
sapply(data, function(x) head(x, 2))
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

# Pairwise comparisions function

``` r
pairwise_perm <- function(c,t,p=1999) {
  out <- NULL
  df <- data$terpene %>% 
    mutate(assess_class_grp = paste(assessment, resistance_class, sep = "-")) %>% 
    filter(mass_type == "dw",
           assess_class_grp == c,
           treatment %in% c("Control", t))%>% 
    select(tree_ID, treatment, assessment, block, compound, mass) %>% 
    pivot_wider(names_from = compound, values_from = mass, values_fill = 0)
  X <- data.frame(
    df %>% select(-treatment, -assessment, -block),
    row.names = 1
  )
  Y <- data.frame(
    df %>% select(tree_ID, treatment, block),
    row.names = 1
  )
  set.seed(146)
  result <- 
    data.frame(
      assess_class_grp = c,
      comparison = paste(sort(unique(df$treatment)), collapse = "_"),
      adonis2(X ~ treatment, data = Y, permutations = p, method = "bray", sqrt.dist = TRUE, strata = Y$block)
    )[1, ]
  out <- rbind(out, result)
  return(out)
}


assess_class_grp <- 
  data$terpene %>% 
  filter(mass_type == "dw") %>% 
  mutate(class_assessment = paste(assessment, resistance_class, sep = "-")) %>% 
  pull(class_assessment) %>% 
  unique()
```

# Results

``` r
list(
  lapply(assess_class_grp, pairwise_perm, t = "EMF"),
  lapply(assess_class_grp, pairwise_perm, t = "FFE"),
  lapply(assess_class_grp, pairwise_perm, t = "FFE+EMF")
) %>% 
  bind_rows() %>% 
  mutate(p_val = `Pr..F.`,
         p_val_adj = round(p.adjust(p_val, "BH"), 4)) %>% 
  rownames_to_column(var = "delete") %>% 
  select(-delete, -`Pr..F.`, -Df, -SumOfSqs) %>%
  separate(assess_class_grp, c("assessment", "resistance_class"), sep = "-") %>% 
  arrange(resistance_class, assessment, comparison) %>% 
  kable(format = "pandoc", caption = "Pairwise comparisons of terpene composition done by permutation (n=1999);\np-values corrected by the Benjamini-Hochberg method.")
```

| assessment | resistance_class | comparison      |        R2 |         F |  p_val | p_val_adj |
|:-----------|:-----------------|:----------------|----------:|----------:|-------:|----------:|
| pre_rust   | MGR              | Control_EMF     | 0.0410018 | 0.9406068 | 0.4305 |    0.5054 |
| pre_rust   | MGR              | Control_FFE     | 0.0381303 | 0.8721202 | 0.4660 |    0.5242 |
| pre_rust   | MGR              | Control_FFE+EMF | 0.0459923 | 1.0606105 | 0.3045 |    0.4952 |
| rust_ctrl  | MGR              | Control_EMF     | 0.2125377 | 4.8582376 | 0.0020 |    0.0180 |
| rust_ctrl  | MGR              | Control_FFE     | 0.1013487 | 2.0300160 | 0.0135 |    0.0729 |
| rust_ctrl  | MGR              | Control_FFE+EMF | 0.0944837 | 1.8781614 | 0.0235 |    0.0906 |
| rust_inoc  | MGR              | Control_EMF     | 0.0545830 | 1.0392177 | 0.3385 |    0.4952 |
| rust_inoc  | MGR              | Control_FFE     | 0.0294264 | 0.5154153 | 0.9370 |    0.9370 |
| rust_inoc  | MGR              | Control_FFE+EMF | 0.0679079 | 1.3113962 | 0.1465 |    0.3596 |
| pre_rust   | QDR              | Control_EMF     | 0.0136975 | 0.9721424 | 0.3920 |    0.5040 |
| pre_rust   | QDR              | Control_FFE     | 0.0104603 | 0.7399614 | 0.7235 |    0.7814 |
| pre_rust   | QDR              | Control_FFE+EMF | 0.0216842 | 1.5515344 | 0.0910 |    0.2457 |
| rust_ctrl  | QDR              | Control_EMF     | 0.0439581 | 2.6208154 | 0.0030 |    0.0203 |
| rust_ctrl  | QDR              | Control_FFE     | 0.0378334 | 2.2806179 | 0.0185 |    0.0832 |
| rust_ctrl  | QDR              | Control_FFE+EMF | 0.0181073 | 1.0695938 | 0.3385 |    0.4952 |
| rust_inoc  | QDR              | Control_EMF     | 0.0170730 | 1.0074344 | 0.3790 |    0.5040 |
| rust_inoc  | QDR              | Control_FFE     | 0.0182796 | 1.0985800 | 0.3175 |    0.4952 |
| rust_inoc  | QDR              | Control_FFE+EMF | 0.0108282 | 0.6349106 | 0.8270 |    0.8588 |
| pre_rust   | susceptible      | Control_EMF     | 0.0340174 | 1.6199053 | 0.0875 |    0.2457 |
| pre_rust   | susceptible      | Control_FFE     | 0.0266680 | 1.2603374 | 0.1810 |    0.3863 |
| pre_rust   | susceptible      | Control_FFE+EMF | 0.0228762 | 1.0769408 | 0.3485 |    0.4952 |
| rust_ctrl  | susceptible      | Control_EMF     | 0.1106859 | 4.8540226 | 0.0005 |    0.0135 |
| rust_ctrl  | susceptible      | Control_FFE     | 0.0853472 | 3.5458182 | 0.0020 |    0.0180 |
| rust_ctrl  | susceptible      | Control_FFE+EMF | 0.0492522 | 1.9685403 | 0.0415 |    0.1401 |
| rust_inoc  | susceptible      | Control_EMF     | 0.0282840 | 1.1060748 | 0.2895 |    0.4952 |
| rust_inoc  | susceptible      | Control_FFE     | 0.0248817 | 0.9696325 | 0.4285 |    0.5054 |
| rust_inoc  | susceptible      | Control_FFE+EMF | 0.0337154 | 1.3258888 | 0.1860 |    0.3863 |

Pairwise comparisons of terpene composition done by permutation
(n=1999); p-values corrected by the Benjamini-Hochberg method.

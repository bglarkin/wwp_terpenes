Identifying indicator terpenes for assessments, treatments, and
resistance classes
================
Beau Larkin

Last updated: 08 November, 2023

- [Description](#description)
- [Package and library installation](#package-and-library-installation)
- [Data](#data)
- [Functions](#functions)
  - [Pre-rust function](#pre-rust-function)
  - [Post-rust function](#post-rust-function)
- [Results](#results)
  - [Pre-rust results](#pre-rust-results)
    - [Plot of indicators and confidence
      intervals](#plot-of-indicators-and-confidence-intervals)
  - [Post-rust results](#post-rust-results)
    - [Rust control seedlings](#rust-control-seedlings)
    - [Plot of indicators and confidence
      intervals](#plot-of-indicators-and-confidence-intervals-1)
    - [Rust-inoculated seedlings](#rust-inoculated-seedlings)
    - [Plot of indicators and confidence
      intervals](#plot-of-indicators-and-confidence-intervals-2)
    - [Heatmaps of terpenes and
      indicators](#heatmaps-of-terpenes-and-indicators)
- [References](#references)

# Description

We’d like to know which terpene compounds associate strongly with
particular assessments or treatments. In multivariate analysis jargon,
the terpenes are “species”, and performing an indicator species analysis
is a creative, but appropriate way to test the specificity and fidelity
of terpenes. Indicators are determined using the concepts “specificity”
and “fidelity”, as explained by [Borcard et
al. 2018](https://doi.org/10.1007/978-3-319-71404-2) (page 120):

> “Indicator value indices are based on the concepts of specificity
> (highest when the species is present in the target group but not
> elsewhere) and fidelity (highest when the species is present in all
> sites of the target group). A high indicator value is obtained by a
> combination of high specificity and fidelity.

The package [indicspecies](https://doi.org/10.1890/08-1823.1) (De
Caceres & Legendre 2009) is used to conduct the indicator species
analysis. Indicspecies identifies species with specificity and fidelity
to sites (in this case, seedlings) grouped by treatments, and then pools
groups, looking for indicators of two, then three, or more (if present)
treatment groups. The groupings can sometimes be difficult to interpret;
for example, when indicators are found for groupings of control and
treatment seedlings.

Finally, the function `strassoc()` is used to produce bootstrapped
confidence intervals on indicators’ strength of association to groups.
The additional post-hoc test reduces the need for or concern over
p-value corrections to `multipatt()`. The bootstrapped statistics and
confidence intervals also allow the production of figures using
`ggplot()` for easy interpretation.

Indicator species analysis here is run on subsets of the seedling
response data:

1.  **Pre-rust inoculation seedlings with resistance classes run
    independently.** This test looks for constitutive differences among
    seedlings and early responses to treatment inoculations.
2.  **Post-rust inoculation seedlings, separated into rust_trt and
    rust_ctrl groups with resistance classes run independently.** This
    test looks for constitutive differences among seedlings (rust_ctrl),
    later responses to treatment inoculations, and induced responses
    crossed with treatments in rust_inoc.
3.  **Post rust inoculated seedlings, treatment controls only, with
    resistance classes run independently.** This test looks for induced
    differences among resistance classes due only to blister rust
    inoculation.

# Package and library installation

Note that messages and code are often hidden in this notebook for
brevity.

``` r
# Package and library installation
packages_needed <- c("tidyverse", "knitr", "indicspecies", "colorspace")
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

# Functions

Two wrapper functions are used to produce summaries of the functions
`multipatt()` and `strassoc()`. Two functions are needed because with
the pre-rust assessment, only treatments are considered within each
resistance class. Post-rust, assessments and treatments must be
considered within each resistance class.

For the permutation test in `multipatt()`, 1999 permutations are used.
In the strength of association test, 1999 bootstrap replicates are
selected in `strassoc()`. These can be changed using the appropriate
arguments in the following functions.

## Pre-rust function

``` r
indVal_prerust_ci <- data.frame()
indVal_pre_pvals <- data.frame()
indic_pre <- function(rc, a="pre_rust", p=1999, nb=1999) {
  df <- data$terpene %>%
    filter(mass_type == "dw",
           assessment == a,
           resistance_class == rc) %>%
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
  indVal <- multipatt(
    X,
    Y$treatment,
    control = how(nperm = p)
  )
  
  ind_compounds <- indVal$sign %>% 
    filter(!is.na(p.value)) %>% 
    rownames()
  
  indVal_boot <- strassoc(X, Y$treatment, func = "IndVal.g", nboot=nb)
  
  indVal_prerust_ci <<- 
    rbind(
      indVal_prerust_ci,
      lapply(indVal_boot, function(x) {
        data.frame(x) %>% 
          rownames_to_column(var = "compound") %>% 
          filter(compound %in% ind_compounds) %>% 
          mutate(resistance_class = rc,
                 assessment = a) %>% 
          select(resistance_class, assessment, compound, everything())
      }) %>% 
        bind_rows(.id = "parameter") %>% 
        pivot_longer(cols = Control:FFE.EMF, names_to = "treatment") %>% 
        mutate(treatment = case_match(treatment, "EMF" ~ "SUIL", "FFE" ~ "META", "FFE.EMF" ~ "MIX", .default = treatment)) %>% 
        pivot_wider(names_from = parameter, values_from = value)
    )
  
  indVal_pre_pvals <<-
    rbind(indVal_pre_pvals, 
          indVal$sign %>% 
            rownames_to_column(var = "compound") %>% 
            mutate(resistance_class = rc,
                   assessment = a,
                   p_val = case_match(p.value, NA ~ 1, .default = p.value),
                   corr_p_val = p.adjust(p_val, method = "BH")) %>% 
            filter(compound %in% ind_compounds, p.value <= 0.05) %>% 
            pivot_longer(cols = 2:5, names_to = "treatment", values_to = "present") %>% 
            filter(present == 1, corr_p_val <= 0.05) %>% 
            mutate(treatment = case_match(treatment, "s.EMF" ~ "SUIL", "s.FFE" ~ "META", "s.FFE+EMF" ~ "MIX", "s.Control" ~ "Control")) %>% 
            select(-index, -present)
    )
  
  print(rc)
  summary(indVal, indvalcomp = TRUE)
  
}
```

## Post-rust function

``` r
indVal_postrust_ci <- data.frame()
indVal_pvals <- data.frame()
indic_post <- function(rc, a, p=1999, nb=1999) {
  df <- data$terpene %>%
    filter(mass_type == "dw",
           assessment == a,
           resistance_class == rc) %>% 
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
  indVal <- multipatt(
    X,
    Y$treatment,
    control = how(nperm = p)
  )
  
  ind_compounds <- indVal$sign %>% 
    filter(!is.na(p.value)) %>% 
    rownames()
  
  indVal_boot <- strassoc(X, Y$treatment, func = "IndVal.g", nboot=nb)
  
  indVal_postrust_ci <<- 
    rbind(
      indVal_postrust_ci,
      lapply(indVal_boot, function(x) {
        data.frame(x) %>% 
          rownames_to_column(var = "compound") %>% 
          filter(compound %in% ind_compounds) %>% 
          mutate(resistance_class = rc,
                 assessment = a) %>% 
          select(resistance_class, assessment, compound, everything())
      }) %>% 
        bind_rows(.id = "parameter") %>% 
        pivot_longer(cols = Control:FFE.EMF, names_to = "treatment") %>% 
        mutate(treatment = case_match(treatment, "EMF" ~ "SUIL", "FFE" ~ "META", "FFE.EMF" ~ "MIX", .default = treatment)) %>% 
        pivot_wider(names_from = parameter, values_from = value)
    )
  
  indVal_pvals <<-
    rbind(indVal_pvals, 
          indVal$sign %>% 
            rownames_to_column(var = "compound") %>% 
            mutate(resistance_class = rc,
                   assessment = a,
                   p_val = case_match(p.value, NA ~ 1, .default = p.value),
                   corr_p_val = p.adjust(p_val, method = "BH")) %>% 
            filter(compound %in% ind_compounds, p.value <= 0.05) %>% 
            pivot_longer(cols = 2:5, names_to = "treatment", values_to = "present") %>% 
            filter(present == 1, corr_p_val <= 0.05) %>% 
            mutate(treatment = case_match(treatment, "s.EMF" ~ "SUIL", "s.FFE" ~ "META", "s.FFE+EMF" ~ "MIX", "s.Control" ~ "Control")) %>% 
            select(-index, -present)
    )
  
  print(paste(a, rc, sep = ", "))
  summary(indVal, indvalcomp = TRUE)
  
}
```

# Results

## Pre-rust results

#### Indicators in QDR seedlings

``` r
indic_pre("QDR")
```

    ## [1] "QDR"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 23
    ##  Selected number of species: 0 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 0 
    ##  Number of species associated to 3 groups: 0 
    ## 
    ##  List of species associated to each combination: 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No terpenes identified.

#### Indicators in susceptible seedlings

``` r
indic_pre("susceptible")
```

    ## [1] "susceptible"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 23
    ##  Selected number of species: 1 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 0 
    ##  Number of species associated to 3 groups: 1 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group EMF+FFE+FFE+EMF  #sps.  1 
    ##              A      B  stat p.value    
    ## abietic 0.8907 1.0000 0.944   5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Abietic acid identified as an indicator in all treatment seedlings, as
opposed to controls, with decent specificity and very high fidelity.

#### Indicators in MGR seedlings

``` r
indic_pre("MGR")
```

    ## [1] "MGR"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 23
    ##  Selected number of species: 1 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 0 
    ##  Number of species associated to 3 groups: 1 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group EMF+FFE+FFE+EMF  #sps.  1 
    ##              A      B  stat p.value    
    ## abietic 0.9073 0.9583 0.932   5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Abietic acid identified as an indicator in all treatment seedlings, as
opposed to controls, with decent specificity and very high fidelity.

#### Summary

**Abietic acid is a consistent indicator of symbiont-treated seedlings
across the board in the pre_rust assessment.**

### Plot of indicators and confidence intervals

The plot below shows indicator statistics and confidence intervals on
single-group comparisons. The statistic shown may not match a
significant pooled-group statistic if one was found using `multipatt()`.
Confidence intervals are based on boostrap replication in `strassoc()`
(n=1999). Confidence intervals which overlap zero mean that the
statistic is non-significant.

![](terpenes_indicators_files/figure-gfm/indVal_prerust_plot-1.png)<!-- -->

## Post-rust results

### Rust control seedlings

Indicator species analyses in rust controls often show terpenes pooled
in groups that combine control and symbiont treatments, which is
difficult to interpret. The first run looked for indicators in single
treatment groups; with none found, the analysis was restricted to
indicators in treatments vs. controls.

#### Indicators in QDR seedlings

``` r
indic_post("QDR", "rust_ctrl")
```

    ## [1] "rust_ctrl, QDR"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 26
    ##  Selected number of species: 7 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 4 
    ##  Number of species associated to 3 groups: 3 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group Control+FFE+EMF  #sps.  2 
    ##                  A      B  stat p.value    
    ## ocimene     0.9900 0.9667 0.978   5e-04 ***
    ## a_terpineol 0.8580 0.9667 0.911   5e-04 ***
    ## 
    ##  Group EMF+FFE+EMF  #sps.  2 
    ##                A      B  stat p.value    
    ## abietic   0.9366 0.9322 0.934   5e-04 ***
    ## palustric 0.8937 0.9153 0.904   5e-04 ***
    ## 
    ##  Group EMF+FFE+FFE+EMF  #sps.  3 
    ##                     A      B  stat p.value    
    ## neoabietic     0.9759 0.9326 0.954   5e-04 ***
    ## levopiramic    0.9724 0.8989 0.935   5e-04 ***
    ## dehydroabietic 1.0000 0.7978 0.893   5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No terpenes are indicators for one treatment group. Neoabietic,
levopiramic, and dehydroabietic are indicators for treatments, as
opposed to controls.

#### Indicators in susceptible seedlings

``` r
indic_post("susceptible", "rust_ctrl")
```

    ## [1] "rust_ctrl, susceptible"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 26
    ##  Selected number of species: 7 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 3 
    ##  Number of species associated to 3 groups: 4 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group Control+FFE+EMF  #sps.  2 
    ##                  A      B  stat p.value    
    ## ocimene     0.9818 1.0000 0.991   5e-04 ***
    ## a_terpineol 0.8562 1.0000 0.925   5e-04 ***
    ## 
    ##  Group EMF+FFE+EMF  #sps.  1 
    ##              A      B  stat p.value    
    ## abietic 0.9136 0.9024 0.908   5e-04 ***
    ## 
    ##  Group EMF+FFE+FFE+EMF  #sps.  4 
    ##                     A      B  stat p.value    
    ## levopiramic    0.9686 0.9180 0.943   5e-04 ***
    ## neoabietic     0.9478 0.8852 0.916   5e-04 ***
    ## palustric      0.9442 0.7869 0.862   5e-04 ***
    ## dehydroabietic 1.0000 0.7377 0.859   5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No terpenes are indicators for one treatment group. Levopiramic,
neoabietic, palustric, and dehydroabietic are indicators for treatments,
as opposed to controls.

#### Indicators in MGR seedlings

``` r
indic_post("MGR", "rust_ctrl")
```

    ## [1] "rust_ctrl, MGR"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 26
    ##  Selected number of species: 6 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 4 
    ##  Number of species associated to 3 groups: 2 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group Control+FFE+EMF  #sps.  2 
    ##                  A      B  stat p.value    
    ## ocimene     0.9735 1.0000 0.987   5e-04 ***
    ## a_terpineol 0.8177 1.0000 0.904   1e-03 ***
    ## 
    ##  Group EMF+FFE+EMF  #sps.  2 
    ##                A      B  stat p.value    
    ## palustric 0.9671 1.0000 0.983   5e-04 ***
    ## abietic   0.9470 0.9000 0.923   5e-04 ***
    ## 
    ##  Group EMF+FFE+FFE+EMF  #sps.  2 
    ##                  A      B  stat p.value    
    ## neoabietic  0.9907 0.9667 0.979   5e-04 ***
    ## levopiramic 0.9807 0.9333 0.957   5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No terpenes are indicators for one treatment group. Levopiramic and
neoabietic are indicators for treatments, as opposed to controls.

#### Summary

- **Constitutive terpenes after rust inoculation segregate along whether
  seedlings were treated with symbionts. There is minimal effect of
  resistance class. Levopiramic, neoabietic, and dehydroabietic are
  consistent terpenes in this group.**
- **Abietic and palustric acids were associated with EMF (both EMF and
  EMF+FFE treatments) across resistance classes among the rust
  controls.**

### Plot of indicators and confidence intervals

The plot below shows indicator statistics and confidence intervals on
single-group comparisons. The statistic shown may not match a
significant pooled-group statistic if one was found using `multipatt()`.
Confidence intervals are based on boostrap replication in `strassoc()`
(n=1999). Confidence intervals which overlap zero mean that the
statistic is non-significant. Zero-overlapping CIs are shown as red on
the plot.

![](terpenes_indicators_files/figure-gfm/indVal_rustctrl_plot-1.png)<!-- -->

### Rust-inoculated seedlings

#### Indicators in QDR seedlings

``` r
indic_post("QDR", "rust_inoc")
```

    ## [1] "rust_inoc, QDR"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 26
    ##  Selected number of species: 3 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 0 
    ##  Number of species associated to 3 groups: 3 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group EMF+FFE+FFE+EMF  #sps.  3 
    ##                  A      B  stat p.value    
    ## ocimene     0.9195 0.9780 0.948   5e-04 ***
    ## a_terpineol 0.9054 0.9560 0.930   5e-04 ***
    ## abietic     0.9941 0.8022 0.893   5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No terpenes are indicators for one treatment group. Ocimene,
a_terpineol, and abietic are associated with all treatments, as opposed
to controls.

#### Indicators in susceptible seedlings

``` r
indic_post("susceptible", "rust_inoc")
```

    ## [1] "rust_inoc, susceptible"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 26
    ##  Selected number of species: 3 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 0 
    ##  Number of species associated to 3 groups: 3 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group EMF+FFE+FFE+EMF  #sps.  3 
    ##                A      B  stat p.value    
    ## ocimene   0.9673 1.0000 0.984   5e-04 ***
    ## palustric 0.9375 0.9833 0.960   5e-04 ***
    ## abietic   0.9805 0.7833 0.876   5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No terpenes are indicators for one treatment group. Ocimene, palustric,
and abietic are associated with all treatments, as opposed to controls.

#### Indicators in MGR seedlings

``` r
indic_post("MGR", "rust_inoc")
```

    ## [1] "rust_inoc, MGR"
    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 26
    ##  Selected number of species: 2 
    ##  Number of species associated to 1 group: 0 
    ##  Number of species associated to 2 groups: 1 
    ##  Number of species associated to 3 groups: 1 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group EMF+FFE  #sps.  1 
    ##              A      B  stat p.value   
    ## abietic 0.8550 0.8947 0.875  0.0075 **
    ## 
    ##  Group EMF+FFE+FFE+EMF  #sps.  1 
    ##             A     B  stat p.value    
    ## ocimene 0.977 1.000 0.988   5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No terpenes are indicators for one treatment group. Only ocimene is
associated with symbiont treatments, as opposed to controls.

#### Summary

- **Terpenes in the post-rust, induced seedlings also often segregate
  along symbiont treatments vs. controls. Consistent indicators in this
  group include ocimene and abietic acid.**
- **Results among resistance classes were spottier here. Where they were
  identified, indicators’ patterns were similar among resistance
  classes, but in several cases, indicators were only identified for one
  or two resistance classes.**
- **Note that rust inoculation is a big hammer on terpenes. An indicator
  analysis performed on treatment=control and rust_ctrl vs. rust_inoc
  assessments revealed that 21 of 26 terpenes were indicators, most of
  rust_inoc (not shown). This didn’t seem an interesting result given
  the lethality of this disease.**

### Plot of indicators and confidence intervals

The plot below shows indicator statistics and confidence intervals on
single-group comparisons. The statistic shown may not match a
significant pooled-group statistic if one was found using `multipatt()`.
Confidence intervals are based on boostrap replication in `strassoc()`
(n=1999). Confidence intervals which overlap zero mean that the
statistic is non-significant. Zero-overlapping CIs are shown as red on
the plot.

![](terpenes_indicators_files/figure-gfm/indVal_rustinoc_plot-1.png)<!-- -->

### Heatmaps of terpenes and indicators

#### Post-rust assessment

The following figure shows variation in terpene masses across
treatments, resistance classes, and post-rust assessments. To create the
heatmap, terpene masses for technical replicates were averaged within
each combination of experimental factors. Then, these averages were
log-transformed to spread the distributions and improve the visual
interpretation of the color gradient. Terpene masses vary among terpene
classes, making low or high molecular weight classes squash to one side
of the color gradient and making differences between cells hard to see.
To take advantage of a full color gradient for each terpene class,
Log-transformed averages were then scaled within each terpene class by
dividing each value by the maximum in that class. This results in values
ranging from 0-1 in each terpene class, displaying the full color
gradient and improving visual interpretation.

Boxes around cells in the heatmap outline indicate that the terpene is a
significant indicator for that combination of experimental factors at
p\<0.05. Indicator statistics were determined using `multipatt()` from
package [indicspecies](http://sites.google.com/site/miqueldecaceres/)
(De Caceres & Legendre 2009) with 1999 permutations. P values were
corrected for multiple comparisons. Significance was visually
corroborated using `strassoc()`, also from package indicspecies, where
95% confidence intervals around the indicator statistic were computed
with 1999 bootstrap samples. Significance is inferred when the 95%
confidence intervals do not include zero.

Note: the indicator statistic returned by `multipatt()` is based on
grouping treatment classes. It does not equal the indicator statistic
returned by `strassoc()`. Either could be used here. Generally, the
statistic returned by `multipatt()` is more robust due to the
availability of treatment groupings, but the confidence intervals
produced by `strassoc()` aren’t limited by the lack of p value
correction and produce a better visual display of differences among
treatments.

**Data wrangling for heatmaps**

Data wrangling is shown here because much of the source data is modified
to produce the figure.

Graphics outputs (not shown) are doubled in size and have all object
sizes doubled as well to improve clarity when shrunk, as per publisher
guidelines.

``` r
terpene_heatmap_data <- 
  data$terpene %>%
  filter(mass_type == "dw", assessment != "pre_rust") %>%
  mutate(
    treatment = case_match(treatment, "EMF" ~ "SUIL", "FFE" ~ "META", "FFE+EMF" ~ "MIX", .default = treatment)) %>% 
  group_by(treatment, assessment, resistance_class, class, compound) %>% 
  summarize(mass = log1p(mean(mass)), .groups = "drop") %>%
  # Terpene masses are averaged and log-transformed in cells to improve visual interpretation.
  left_join(
    indVal_pvals %>%
      mutate(sig = 0.5),
    # A continuous placeholder variable must be created so ggplot2 can draw significance boxes later.
    by = join_by(treatment, assessment, resistance_class, compound)
  ) %>% 
  mutate(
    assessment = case_match(assessment, "rust_ctrl" ~ "Pathogen-", "rust_inoc" ~ "Pathogen+", .default = assessment),
    resistance_class = case_match(resistance_class, "susceptible" ~ "Susceptible", .default = resistance_class),
    resistance_class = factor(resistance_class, ordered = TRUE, levels = c("QDR", "Susceptible", "MGR")),
    treatment = factor(treatment, levels = c("Control", "SUIL", "META", "MIX"), ordered = TRUE),
    class = case_match(class, "diterpene" ~ "Diterpene acids", "monoterpene" ~ "Monoterpene", "sesquiterpene" ~ "Sesquiterpene"),
    compound = factor(compound)
  ) %>% 
  group_by(class) %>%
  mutate(mass_scl = mass/max(mass)) %>%
  # Raw terpene masses vary among terpene classes, making differences between cells hard to see. 
  # Log-transformed averages are scaled within each terpene class to improve visual representation.
  ungroup()
```

![](terpenes_indicators_files/figure-gfm/indVal_heatmap_plot-1.png)<!-- -->

**Figure, Sampling Period 2:** Heatmap of 26 terpenes identified in
foliar tissue of P. monticola with quantitative disease resistance,
susceptibility to disease, or major gene resistance inoculated with
fungal symbionts. Vertical panels represent Pathogen- and Pathogen+
treatments in all resistance classes at collection two. Horizontal
panels are separated according to terpene class. Shading in each cell
relates to terpene concentration, which is an average of individual
seedling values in each cell (n = 9-30) . Averages were log-transformed
improve contrast and enhance visual interpretation of the figure.
Log-transformed values were then divided by the maximum value of each
terpene class to relativize the color gradient across terpene classes,
again to enhance visual interpretation of the figure. Analysis and
significance tests were done on untransformed data as described in the
methods. Cells with solid black outlines indicate significant
associations at p \< 0.05 (values corrected for multiple tests) as
determined by indicator species analysis.

#### Pre-rust assessment

Included for a supplemental figure

``` r
terpene_pre_heatmap_data <- 
  data$terpene %>%
  filter(mass_type == "dw", assessment == "pre_rust") %>% 
  mutate(
    treatment = case_match(treatment, "EMF" ~ "SUIL", "FFE" ~ "META", "FFE+EMF" ~ "MIX", .default = treatment)) %>% 
  group_by(treatment, resistance_class, class, compound) %>% 
  summarize(mass = log1p(mean(mass)), .groups = "drop") %>%
  # Terpene masses are averaged and log-transformed in cells to improve visual interpretation.
  left_join(
    indVal_pre_pvals %>%
      mutate(sig = 0.5),
    # A continuous placeholder variable must be created so ggplot2 can draw significance boxes later.
    by = join_by(treatment, resistance_class, compound)
  ) %>% 
  mutate(
    resistance_class = case_match(resistance_class, "susceptible" ~ "Susceptible", .default = resistance_class),
    resistance_class = factor(resistance_class, ordered = TRUE, levels = c("QDR", "Susceptible", "MGR")),
    treatment = factor(treatment, levels = c("Control", "SUIL", "META", "MIX"), ordered = TRUE),
    class = case_match(class, "diterpene" ~ "Diterpene acids", "monoterpene" ~ "Monoterpene", "sesquiterpene" ~ "Sesquiterpene"),
    compound = factor(compound)
  ) %>% 
  group_by(class) %>%
  mutate(mass_scl = mass/max(mass)) %>%
  # Raw terpene masses vary among terpene classes, making differences between cells hard to see. 
  # Log-transformed averages are scaled within each terpene class to improve visual representation.
  ungroup()
```

![](terpenes_indicators_files/figure-gfm/indVal_pre_heatmap_plot-1.png)<!-- -->

**Figure, Sampling Period 1:** Heatmap of 23 terpenes identified in
foliar tissue of P. monticola with quantitative disease resistance,
susceptibility to disease, or major gene resistance inoculated with
fungal symbionts. Vertical panels represent Pathogen- and Pathogen+
treatments in all resistance classes at sampling period 2. Horizontal
panels are separated according to terpene class. Shading in each cell
relates to terpene concentration, which is an average of individual
seedling values in each cell (n = 8-48) . Averages were log-transformed
improve contrast and enhance visual interpretation of the figure.
Log-transformed values were then divided by the maximum value of each
terpene class to relativize the color gradient across terpene classes,
again to enhance visual interpretation of the figure. Analysis and
significance tests were done on untransformed data as described in the
methods. Cells with solid black outlines indicate significant
associations at p \< 0.05 (values corrected for multiple tests) as
determined by indicator species analysis.

# References

``` r
print(citation("indicspecies"), bibtex = FALSE)
```

    ## 
    ## To cite 'indicspecies' package in publications use:
    ## 
    ##   De Caceres, M., Legendre, P. (2009). Associations between species and
    ##   groups of sites: indices and statistical inference. Ecology, URL
    ##   http://sites.google.com/site/miqueldecaceres/
    ## 
    ## Thank you for using 'indicspecies'

``` r
print(citation("tidyverse"), bibtex = FALSE)
```

    ## 
    ## To cite package 'tidyverse' in publications use:
    ## 
    ##   Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R,
    ##   Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller
    ##   E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V,
    ##   Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). "Welcome to
    ##   the tidyverse." _Journal of Open Source Software_, *4*(43), 1686.
    ##   doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.

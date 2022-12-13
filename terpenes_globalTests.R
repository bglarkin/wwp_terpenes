#' ---
#' title: "Global tests of terpene composition"
#' author: "Beau Larkin"
#' date: "2022-12-08"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#'     fig_width: 7
#'     fig_height: 6
#' ---
#'
#' # Description
#' Permutation tests are used to test multivariate differences in terpene composition
#' among control, treatment, and assessment groups. One test is done for each resistance class.
#' Then, a global test is done using all terpenes and testing resistance classes, treatments,
#' and assessments in a three-way anova by permutation.
#'
#' # Package and library installation
#' Note that messages and code are often hidden in this notebook for brevity.
# Package and library installation
packages_needed <-
  c("tidyverse", "magrittr", "knitr", "vegan", "colorspace")
packages_installed <-
  packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed))
  install.packages(packages_needed[!packages_installed])
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
  library(packages_needed[i], character.only = T)
}
#+ ggstyle
# Load ggplot styles and themes from text file
source("gg_style.txt")
#'
#' # Data
#' See [data_etl.md](https://github.com/bglarkin/wwp_terpenes/blob/main/data_etl.md) for more description of
#' the source data. Header views of each data table are presented here.
#' Names only provided here.
#+ data_source,message=FALSE,results=FALSE
source("data_etl.R")
#+ data_headers
sapply(data, function(x)
  head(x, 2))
#'
#' # Functions
#' The following function produces the permutation tests and visual ordination figures for
#' each resistance class (permutations = 1999). The test is run on individual trees; the ordination figures show
#' centroids and standard errors for assessment and treatment groups. Supplemental ordinations of
#' terpene compounds are also shown.
#'
#' ## Terms
#' - permutations = 19999
#' - data standardization: standardize columns (scale x to zero mean and unit variance)
#' - distance metric: euclidean
#' - ordination: PCA
#+ perm_pcoa_function
terpene_pca <-
  function(c,
           bar_wd = 0.008,
           bar_sz = 0.2,
           pt_sz = 3,
           p = 1999) {
    cat("---------------------------------------------------------------------\n")
    cat(paste("Resistance type", c, "selected."))
    cat("\n---------------------------------------------------------------------\n")
    
    # Load styles inside function
    source("gg_style.txt")
    
    df <- data$terpene %>%
      filter(mass_type == "dw",
             resistance_class == c) %>%
      mutate(tree_key = paste(tree_ID, year, sep = "-")) %>%
      select(tree_key, treatment, assessment, compound, mass) %>%
      pivot_wider(
        names_from = compound,
        values_from = mass,
        values_fill = 0
      )
    X <- data.frame(df %>% select(-treatment,-assessment), row.names = 1)
    Y <- data.frame(df %>% select(treatment, assessment))
    set.seed(123)
    perm_test <-
      adonis2(
        scale(X) ~ assessment * treatment,
        data = Y,
        permutations = p,
        method = "euclidean"
      )
    
    terp_pca <- rda(X, scale = TRUE)
    sites <-
      data.frame(scores(terp_pca, "sites", choices = c(1, 2))) %>%
      rownames_to_column(var = "tree_key") %>%
      left_join(df %>% select(tree_key, treatment, assessment), by = "tree_key")
    site_centers <-
      sites %>%
      group_by(assessment, treatment) %>%
      summarize(
        PC1_mean = mean(PC1),
        PC1_se_pos = PC1_mean + (sd(PC1) / sqrt(length(PC1))),
        PC1_se_neg = PC1_mean - (sd(PC1) / sqrt(length(PC1))),
        PC2_mean = mean(PC2),
        PC2_se_pos = PC2_mean + (sd(PC2) / sqrt(length(PC2))),
        PC2_se_neg = PC2_mean - (sd(PC2) / sqrt(length(PC2))),
        .groups = "drop"
      )
    labs_pct <-
      round((terp_pca$CA$eig / sum(terp_pca$CA$eig))[1:2] * 100, 0)
    terp_centers <-
      data.frame(scores(terp_pca, "species", choices = c(1, 2))) %>%
      rownames_to_column(var = "compound")
    
    bar_width_d1 <-
      with(site_centers, max(PC1_se_pos) - min(PC1_se_neg)) * bar_wd
    bar_width_d2 <-
      with(site_centers, max(PC2_se_pos) - min(PC2_se_neg)) * bar_wd
    
    plot_ord <-
      ggplot(site_centers, aes(x = PC1_mean, y = PC2_mean)) +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_errorbar(
        aes(x = PC1_mean, ymin = PC2_se_neg, ymax = PC2_se_pos),
        width = bar_width_d1,
        size = bar_sz
      ) +
      geom_errorbar(
        aes(y = PC2_mean, xmin = PC1_se_neg, xmax = PC1_se_pos),
        width = bar_width_d2,
        size = bar_sz
      ) +
      geom_point(aes(shape = assessment, fill = treatment),
                 size = pt_sz,
                 stroke = bar_sz) +
      labs(
        x = paste0("Component 1, ", labs_pct[1], "% variation explained"),
        y = paste0("Component 2, ", labs_pct[2], "% variation explained"),
        title = paste0("Terpenes in ", c, " families")
      ) +
      scale_shape_manual(name = "Assessment", values = c(21, 22, 24)) +
      scale_fill_discrete_qualitative(name = "Treatment", palette = "Harmonic") +
      guides(fill = guide_legend(override.aes = list(shape = 21)),
             shape = guide_legend(override.aes = list(fill = "gray50"))) +
      theme_bgl
    
    plot_compounds <-
      ggplot(terp_centers, aes(x = PC1, y = PC2)) +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_label(aes(label = compound),
                 family = "serif",
                 size = 8 * 0.36) +
      labs(
        x = paste0("Component 1, ", labs_pct[1], "% variation explained"),
        y = paste0("Component 2, ", labs_pct[2], "% variation explained"),
        title = paste0("Terpenes in ", c, " families")
      ) +
      theme_bgl
    
    out <- list(permutation_test_result = perm_test)
    
    #+ figure_ordination
    print(plot_ord)
    print(plot_compounds)
    
    return(out)
    
  }
#'
#' # Results
#' ## Tests within resistance classes
#' Use of standardized terpene masses in a euclidean distance matrix greatly improves clustering
#' and increases separation in a permutation test.
#'
#' ## Quantitative gene resistance class seedlings
#+ qdr_test,echo=FALSE
terpene_pca("QDR")
#' ## Susceptible resistance class seedlings
#+ susc_test,echo=FALSE
terpene_pca("susceptible")
#' ## Major gene resistance class seedlings
#+ mgr_test,echo=FALSE
terpene_pca("MGR")
#'
#' ## Global test
#+ global_data
g_df <- data$terpene %>%
  filter(mass_type == "dw") %>%
  mutate(tree_key = paste(tree_ID, year, sep = "-")) %>%
  select(tree_key,
         resistance_class,
         treatment,
         assessment,
         compound,
         mass) %>%
  pivot_wider(names_from = compound,
              values_from = mass,
              values_fill = 0)
g_X <-
  data.frame(g_df %>% select(-resistance_class, -treatment, -assessment),
             row.names = 1)
g_Y <-
  data.frame(g_df %>% select(resistance_class, treatment, assessment))
#+ global_perm_test
set.seed(123)
g_perm_test <-
  adonis2(
    scale(g_X) ~ resistance_class * assessment * treatment,
    data = g_Y,
    permutations = 1999,
    method = "euclidean"
  )
#+ global_perm_test_table,echo=FALSE
g_perm_test

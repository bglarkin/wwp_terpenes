#' ---
#' title: "Global tests of terpene composition"
#' author: "Beau Larkin"
#' date: "2022-12-08"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#'     fig_width: 9.5
#'     fig_height: 9
#' ---
#' 
#' # Description
#' Permutation tests are used to test multivariate differences in terpene composition
#' among control, treatment, and assessment groups. One test is done for each resistance class.
#' 
#' # Package and library installation
#' Note that messages and code are often hidden in this notebook for brevity.
# Package and library installation
packages_needed <- c("tidyverse", "magrittr", "knitr", "vegan", "colorspace")
packages_installed <- packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(! packages_installed))
  install.packages(packages_needed[! packages_installed])
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
  library(packages_needed[i], character.only = T)
}
# Load ggplot styles and themes from text file
#+ ggstyle
source("gg_style.txt")
#'
#' # Data
#' See [data_etl.md](https://github.com/bglarkin/wwp_terpenes/blob/main/data_etl.md) for a brief data dictionary. 
#' Names only provided here.
#+ data_source,message=FALSE,results=FALSE
source("data_etl.R")
#+ data_names
names(data)
#'
#' # Results
#' The following function produces the permutation tests and visual ordination figures for
#' each resistance class (permutations = 1999). The test is run on individual trees; the ordination figures show
#' centroids and standard errors for assessment and treatment groups.
#+ perm_pcoa_function
terpene_pcoa <- function(c, dim1_exp = 1, dim2_exp = 1, bar_wd = 0.008, bar_sz = 0.2, pt_sz = 3) {

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
    pivot_wider(names_from = compound, values_from = mass, values_fill = 0)
  terp <-
    data.frame(df %>% select(-treatment, -assessment), row.names = 1)
  expl <- data.frame(df %>% select(treatment, assessment))
  set.seed(123)
  perm_test <-
    adonis2(
      terp ~ treatment * assessment,
      data = expl,
      permutations = 1999,
      method = "bray",
      sqrt.dist = TRUE,
      by = "terms"
    )

  terp_bray <- vegdist(terp, "bray")
  terp_pcoa <-
    cmdscale(sqrt(terp_bray), k = (nrow(terp) - 1), eig = TRUE)
  sites <-
    data.frame(scores(terp_pcoa, "sites", choices = c(1, 2))) %>%
    rownames_to_column(var = "tree_key") %>%
    left_join(df %>% select(tree_key, treatment, assessment), by = "tree_key")
  site_centers <-
    sites %>%
    group_by(assessment, treatment) %>%
    summarize(Dim1_mean = mean(Dim1),
              Dim1_se_pos = Dim1_mean + (sd(Dim1) / sqrt(length(Dim1))),
              Dim1_se_neg = Dim1_mean - (sd(Dim1) / sqrt(length(Dim1))),
              Dim2_mean = mean(Dim2),
              Dim2_se_pos = Dim2_mean + (sd(Dim2) / sqrt(length(Dim2))),
              Dim2_se_neg = Dim2_mean - (sd(Dim2) / sqrt(length(Dim2))),
              .groups = "drop")
  labs_pct <-
    round((terp_pcoa$eig / sum(terp_pcoa$eig))[1:2] * 100, 0)
  terp_wa <-
    data.frame(wascores(terp_pcoa$points[, 1:2], terp, expand = FALSE)) %>%
    rownames_to_column(var = "compound")

  bar_width_d1 <- with(site_centers, max(Dim1_se_pos) - min(Dim1_se_neg)) * bar_wd
  bar_width_d2 <- with(site_centers, max(Dim2_se_pos) - min(Dim2_se_neg)) * bar_wd

  plot_ord <-
    ggplot(site_centers, aes(x = Dim1_mean, y = Dim2_mean)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(
      aes(x = Dim1_mean, ymin = Dim2_se_neg, ymax = Dim2_se_pos),
      width = bar_width_d1,
      size = bar_sz
    ) +
    geom_errorbar(
      aes(y = Dim2_mean, xmin = Dim1_se_neg, xmax = Dim1_se_pos),
      width = bar_width_d2,
      size = bar_sz
    ) +
    geom_point(
      aes(shape = assessment, fill = treatment),
      size = pt_sz,
      stroke = bar_sz) +
    geom_text(
      data = terp_wa,
      aes(x = X1, y = X2, label = compound),
      family = "serif",
      size = 8 * 0.36
    ) +
    labs(
      x = paste0("Dimension 1, ", labs_pct[1], "% variation explained"),
      y = paste0("Dimension 2, ", labs_pct[2], "% variation explained"),
      title = paste0("Terpenes in ", c, " families")
    ) +
    scale_shape_manual(name = "Assessment", values = c(21,22,24)) +
    scale_fill_discrete_qualitative(name = "Treatment", palette = "Harmonic") +
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = "gray50"))) +
    theme_bgl

  out <- list(
    permutation_test_result = perm_test
  )
  
  #+ figure_ordination
  print(plot_ord)

  return(out)

}
#'
#' ## Susceptible resistance class seedlings
#+ susc_test
terpene_pcoa("susceptible")
#' ## Major gene resistance class seedlings
#+ mgr_test
terpene_pcoa("MGR")
#' ## Quantitative gene resistance class seedlings
#+ qdr_test
terpene_pcoa("QDR")



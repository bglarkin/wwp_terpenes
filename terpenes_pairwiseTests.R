#' ---
#' title: "Pairwise tests of terpene composition"
#' author: "Beau Larkin"
#' date: "2022-12-09"
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
#' Permutation tests are used to test pairwise multivariate differences in terpene composition
#' within resistance classes and assessments, and between treatments. The tests are 
#' permuted within experimental greenhouse blocks. [Benjamini-Hochberg](doi:10.1111/j.2517-6161.1995.tb02031.x) 
#' corrections are used to adjust p values for multiple tests. 
#' 
#' # Package and library installation
#' Note that messages and code are often hidden in this notebook for brevity.
# Package and library installation
packages_needed <- c("tidyverse", "knitr", "vegan", "colorspace")
packages_installed <- packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(! packages_installed))
  install.packages(packages_needed[! packages_installed])
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
sapply(data, function(x) head(x, 2))
#'
#' # Pairwise comparisions function



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
      adonis2(scale(X) ~ treatment, data = Y, permutations = p, method = "euclidean", sqrt.dist = FALSE, strata = Y$block)
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

#' # Results
#' 
#' 
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




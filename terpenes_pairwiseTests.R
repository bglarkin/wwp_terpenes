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
#' # Pairwise comparisons function
#' This function produces pairwise group comparisons between fungal controls and treatments
#' for each pairwise comparison among assessments and resistance classes. A permutation test is 
#' used in each comparison, using a blocks as strata in the model, 
#' and a correction is made to p-values due to the number of comparisons 
#' being done. 
#' 
#' 
#' A dummy variable combining assessments and resistance_classes is needed to reduce the number of nested calls to `lapply()`. 
#' The variable is `assess_class_grp`, below.
#' 
#' ## Terms
#' - permutations = 19999
#' - data standardization: standardize columns (scale x to zero mean and unit variance)
#' - distance metric: euclidean
#+ pairwise_perm_function
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
#' See output table below.
#+ func_apply
func_apply <- 
list(
  lapply(assess_class_grp, pairwise_perm, t = "EMF"),
  lapply(assess_class_grp, pairwise_perm, t = "FFE"),
  lapply(assess_class_grp, pairwise_perm, t = "FFE+EMF")
)
#+ func_result
result_tab <- 
  func_apply %>% 
  bind_rows() %>% 
  mutate(p_val = `Pr..F.`,
         p_val_adj = round(p.adjust(p_val, "BH"), 4),
         sig_05 = case_when(p_val_adj <= 0.05 ~ "*", TRUE ~ "")) %>% 
  rownames_to_column(var = "delete") %>% 
  select(-delete, -`Pr..F.`, -Df, -SumOfSqs) %>%
  separate(assess_class_grp, c("assessment", "resistance_class"), sep = "-") %>% 
  arrange(resistance_class, assessment, comparison)
#+ func_table,echo=FALSE
result_tab %>% 
  kable(format = "pandoc", 
        caption = "Pairwise comparisons of terpene composition done by permutation (n=1999); p-values corrected by the Benjamini-Hochberg method.")
#+ export,echo=FALSE,result=FALSE
write_csv(result_tab, "pairwise_perm_table.csv")



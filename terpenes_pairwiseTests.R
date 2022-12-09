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
#' # Results
#' 
#' 



pairwise_perm <- function(c, t) {
  df <- data$terpene %>% 
    mutate(class_assessment = paste(assessment, resistance_class, sep = "_")) %>% 
    filter(mass_type == "dw",
           class_assessment == c,
           treatment %in% c("Control", t))%>% 
    select(tree_ID, treatment, assessment, block, compound, mass) %>% 
    pivot_wider(names_from = compound, values_from = mass)
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
      resistance_class = c,
      comparison = paste(sort(unique(df$treatment)), collapse = "_"),
      adonis2(X ~ treatment, data = Y, permutations = 99, method = "bray", sqrt.dist = TRUE, strata = Y$block)
    )[1, ]
  out <- rbind(out, result)
  return(out)
}

out <- NULL
pairwise_perm("pre_rust_MGR", "FFE")

# "pre_rust_susceptible"  "pre_rust_MGR"          "pre_rust_QDR"          "rust_inoc_susceptible" "rust_inoc_MGR"         "rust_inoc_QDR"         "rust_ctrl_QDR"        
# "rust_ctrl_MGR"         "rust_ctrl_susceptible"




data$terpene %>% 
  filter(mass_type == "dw") %>% 
  mutate(class_assessment = paste(assessment, resistance_class, sep = "_")) %>% 
  pull(class_assessment) %>% 
  unique()





pairwise_perm_inner <- function(a, c, t) {
  df <- data$terpene %>% 
    filter(mass_type == "dw",
           resistance_class == c, 
           assessment == a,
           treatment %in% c("Control", t)) %>% 
    select(tree_ID, treatment, assessment, block, compound, mass) %>% 
    pivot_wider(names_from = compound, values_from = mass)
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
      resistance_class = c,
      comparison = paste(sort(unique(df$treatment)), collapse = "_"),
      adonis2(X ~ treatment, data = Y, permutations = 99, method = "bray", sqrt.dist = TRUE, strata = Y$block)
    )[1, ]
  out <- rbind(out, result)
  return(out)
}



pairwise_perm_outer <- function(a) {
  
  
  list(
    lapply(classes, pairwise_perm_inner, t = "EMF", a = a),
    lapply(classes, pairwise_perm_inner, t = "FFE", a = a),
    lapply(classes, pairwise_perm_inner, t = "FFE+EMF", a = a)
  ) %>% 
    bind_rows() %>% 
    rownames_to_column(var = "delete") %>% 
    select(-delete, -Df, -SumOfSqs) %>%
    arrange(resistance_class) %>% 
    kable(format = "pandoc")
  
}


out <- NULL
classes <- c("QDR", "susceptible", "MGR")
assessments <- c("pre_rust", "rust_inoc", "rust_ctrl")
pairwise_perm_outer(assessments)

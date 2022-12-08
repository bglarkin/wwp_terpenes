#' ---
#' title: "Extract, transform, and load data"
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
#' # Package and library installation
#' Note that messages and code are often hidden in this notebook for brevity.
# Package and library installation
packages_needed <- c("tidyverse", "magrittr", "knitr")
packages_installed <- packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(! packages_installed))
  install.packages(packages_needed[! packages_installed])
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
  library(packages_needed[i], character.only = T)
}
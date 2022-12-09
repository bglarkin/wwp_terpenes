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
#' # Description
#' This script pulls data from separate files in the **Database** folder, adds
#' a few new variables, and checks data for consistency. This script will be 
#' called from other analysis scripts so that the ETL process isn't repeated, 
#' making the analysis scripts a little shorter.
#' 
#' ## Pre-cleaning
#+ cleaning,message=FALSE
rm(list = ls()) 
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
#'
#' # Data
#' - Data are extracted from the associated directory and loaded into a list.
#' - The `resistance_class` variable is synthesized and added to several tables
#+ data,echo=TRUE,results=TRUE
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
#'
#' ## Data and variable views
#' ### `terpene_meta`
#' Experimental design and metadata on trees that had terpenes extracted. Useful for an **env** explanatory
#' file to go along with a distance matrix produced from `terpene`.
data[1]
#' ### `terpene`
#' Granular terpene masses extracted from a subset of experimental seedlings. Also includes experimental
#' design and metadata for convenience. Normally, only dry weights of terpenes are used in analysis,
#' but wet weights are also included here. 
data[2]
#' ### `tree_height`
#' Only applicable height measurements are included here to simplify use of height data in correlations.
data[3]
#' ### `tree_meta`
#' Experimental design and metadata on all trees (not just those that had terpenes extracted).
data[4]
#' ### `tree_rust_response`
#' Disease response traits on all seedlings.
data[5]

#' ---
#' title: "Identifying indicator terpenes for assessments, treatments, and resistance classes"
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
#' We'd like to know which terpene compounds associate strongly with particular assessments or treatments. 
#' In multivariate analysis jargon, the terpenes are "species", and performing an indicator species
#' analysis is a creative, but appropriate way to test the specificity and fidelity of terpenes. Indicators are determined
#' using the concepts "specificity" and "fidelity", as explained by [Borcard et al. 2018](https://doi.org/10.1007/978-3-319-71404-2) (page 120):
#' 
#' > "Indicator value indices are based on the concepts of specificity (highest when the species is present 
#' > in the target group but not elsewhere) and fidelity (highest when the species is present in all sites of the target group). 
#' > A high indicator value is obtained by a combination of high specificity and fidelity.
#'
#' # Package and library installation
#' Note that messages and code are often hidden in this notebook for brevity.
# Package and library installation
packages_needed <- c("tidyverse", "knitr", "indicspecies")
packages_installed <-
  packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed))
  install.packages(packages_needed[!packages_installed])
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
  library(packages_needed[i], character.only = T)
}
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
#' Two wrapper functions are used to produce summaries of the function `multipatt()`. Two functions
#' are needed because with the pre-rust assessment, only treatments are considered within each 
#' resistance class. Post-rust, assessments and treatments must be considered within each resistance class.
#' ## Pre-rust function 
#+ pre-rust function
indic_pre <- function(rc) {
  df <- data$terpene %>%
    filter(mass_type == "dw",
           assessment == "pre_rust",
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
    control = how(nperm = 999)
  )
  print(rc)
  summary(indVal, indvalcomp = TRUE)
  
}
#'
#' ## Post-rust function
#+ post_rust_function
indic_post <- function(rc, a) {
  out <- NULL
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
    control = how(nperm = 999)
  )
  print(paste(a, rc, sep = ", "))
  summary(indVal, indvalcomp = TRUE)
  
}
#' 
#' # Results
#' ## Pre-rust results
#' #### Indicators in QDR seedlings
#+ qdr_pre
indic_pre("QDR")
#' No terpenes identified. 
#' #### Indicators in susceptible seedlings
#+ sus_pre
indic_pre("susceptible")
#' Abietic acid identified as an indicator in all treatment seedlings, as opposed to controls, 
#' with decent specificity and very high fidelity.
#' #### Indicators in MGR seedlings
#+ mgr_pre
indic_pre("MGR")
#' Abietic acid identified as an indicator in all treatment seedlings, as opposed to controls, 
#' with decent specificity and very high fidelity.
#' 
#' ## Post-rust results
#' ### Rust control seedlings
#' Indicator species analyses in rust controls often show terpenes pooled in groups that combine 
#' control and symbiont treatments, which is difficult to interpret.
#' #### Indicators in QDR seedlings
#+ qdr_post_ctrl
indic_post("QDR", "rust_ctrl")
#' 
#' #### Indicators in susceptible seedlings
#+ sus_post_ctrl
indic_post("susceptible", "rust_ctrl")
#' 
#' #### Indicators in MGR seedlings
#+ mgr_post_ctrl
indic_post("MGR", "rust_ctrl")
#' 
#' ### Rust-inoculated seedlings
#' #### Indicators in QDR seedlings
#+ qdr_post_inoc
indic_post("QDR", "rust_inoc")
#' 
#' #### Indicators in susceptible seedlings
#+ sus_post_inoc
indic_post("susceptible", "rust_inoc")
#' 
#' #### Indicators in MGR seedlings
#+ mgr_post_inoc
indic_post("MGR", "rust_inoc")
#' 
#' 

# Consider a way to look at differences based on rust inoculation...in all resistance classes but only control treatments
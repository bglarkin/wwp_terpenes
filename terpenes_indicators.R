#' ---
#' title: "Identifying indicator terpenes for assessments, treatments, and resistance classes"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
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
#' The package [indicspecies](https://doi.org/10.1890/08-1823.1) (De Caceres & Legendre 2009) is used to 
#' conduct the indicator species analysis. Indicspecies identifies species with specificity and fidelity to sites
#' (in this case, seedlings) grouped by treatments, and then pools groups, looking for indicators of two, then three, 
#' or more (if present) treatment groups. The groupings can 
#' sometimes be difficult to interpret; for example, when indicators are found for groupings of control and treatment
#' seedlings. 
#' 
#' Finally, the function `strassoc()` is used to produce bootstrapped confidence intervals on indicators' strength of 
#' association to groups. The additional post-hoc test hopefully reduces the need for or concern over p-value corrections to 
#' `multipatt()`. The bootstrapped statistics and confidence intervals also allow the production of figures using
#' `ggplot()` for easy interpretation. 
#' 
#' Indicator species analysis here is run on subsets of the seedling response data:
#' 
#' 1. **Pre-rust inoculation seedlings with resistance classes run independently.** This test looks for constitutive differences
#' among seedlings and early responses to treatment inoculations. 
#' 1. **Post-rust inoculation seedlings, separated into rust_trt and rust_ctrl groups with resistance classes run independently.** 
#' This test looks for constitutive differences among seedlings (rust_ctrl), later responses to treatment inoculations, and
#' induced responses crossed with treatments in rust_inoc. 
#' 1. **Post rust inoculated seedlings, treatment controls only, with resistance classes run independently.** This test looks for induced differences among resistance classes 
#' due only to blister rust inoculation. 
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
#' Two wrapper functions are used to produce summaries of the functions `multipatt()` and `strassoc()`. Two functions
#' are needed because with the pre-rust assessment, only treatments are considered within each 
#' resistance class. Post-rust, assessments and treatments must be considered within each resistance class.
#' 
#' ## Pre-rust function 
#+ pre-rust function
indVal_prerust_ci <- data.frame()
indic_pre <- function(rc, a="pre_rust", p=999, nb=999) {
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
    pivot_wider(names_from = parameter, values_from = value)
  )
  
  print(rc)
  summary(indVal, indvalcomp = TRUE)
  
}
#'
#' ## Post-rust function
#+ post_rust_function
indVal_postrust_ci <- data.frame()
indic_post <- function(rc, a, p=999, nb=999) {
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
        pivot_wider(names_from = parameter, values_from = value)
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
#' 
#' #### Indicators in susceptible seedlings
#+ sus_pre
indic_pre("susceptible")
#' Abietic acid identified as an indicator in all treatment seedlings, as opposed to controls, 
#' with decent specificity and very high fidelity.
#' 
#' #### Indicators in MGR seedlings
#+ mgr_pre
indic_pre("MGR")
#' Abietic acid identified as an indicator in all treatment seedlings, as opposed to controls, 
#' with decent specificity and very high fidelity.
#' 
#' #### Summary
#' **Abietic acid is a consistent indicator of symbiont-treated seedlings across the board
#' in the pre_rust assessment.**
#' 
#' ### Plot of indicators and confidence intervals
#' The plot below shows indicator statistics and confidence intervals on single-group comparisons. 
#' The statistic shown may not match a significant pooled-group statistic if one was found using 
#' `multipatt()`. Confidence intervals are based on boostrap replication in `strassoc()` (n=1000). 
#' Confidence intervals which overlap zero mean that the statistic is non-significant. 
#+ ggstyle,echo=FALSE
source("gg_style.txt")
#+ indVal_prerust_plot,echo=FALSE,fig.dim=c(9,6)
indVal_prerust_ci %>%
ggplot(aes(x = treatment, y = stat)) +
  facet_grid(rows = vars(compound), cols = vars(resistance_class)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI)) +
  labs(x = "", 
       y = "Single-group indicator statistic", 
       title = "Indicator terpenes in pre-rust inoculation seedlings") +
  theme_bw() +
  theme_bgl
#' 
#' ## Post-rust results
#' 
#' ### Rust control seedlings
#' Indicator species analyses in rust controls often show terpenes pooled in groups that combine 
#' control and symbiont treatments, which is difficult to interpret. The first run looked for indicators in
#' single treatment groups; with none found, the analysis was restricted to indicators in treatments vs. controls.
#' 
#' #### Indicators in QDR seedlings
#+ qdr_post_ctrl
indic_post("QDR", "rust_ctrl")
#' No terpenes are indicators for one treatment group. 
#' Neoabietic, levopiramic, and dehydroabietic are indicators for treatments, as opposed to controls.
#' 
#' #### Indicators in susceptible seedlings
#+ sus_post_ctrl
indic_post("susceptible", "rust_ctrl")
#' No terpenes are indicators for one treatment group. 
#' Levopiramic, neoabietic, palustric, and dehydroabietic are indicators for treatments, as opposed to controls. 
#' 
#' #### Indicators in MGR seedlings
#+ mgr_post_ctrl
indic_post("MGR", "rust_ctrl")
#' No terpenes are indicators for one treatment group. 
#' Levopiramic and neoabietic are indicators for treatments, as opposed to controls. 
#' 
#' #### Summary
#' - **Constitutive terpenes after rust inoculation segregate along whether seedlings were treated with symbionts. There is minimal effect
#' of resistance class. Levopiramic, neoabietic, and dehydroabietic are consistent terpenes in this group.**
#' - **Abietic and palustric acids were associated with EMF (both EMF and EMF+FFE treatments) across resistance classes
#' among the rust controls.**
#' 
#' ### Plot of indicators and confidence intervals
#' The plot below shows indicator statistics and confidence intervals on single-group comparisons. 
#' The statistic shown may not match a significant pooled-group statistic if one was found using 
#' `multipatt()`. Confidence intervals are based on boostrap replication in `strassoc()` (n=1000). 
#' Confidence intervals which overlap zero mean that the statistic is non-significant. Zero-overlapping
#' CIs are shown as red on the plot. 
#+ indVal_rustctrl_plot,echo=FALSE,fig.dim=c(9,16)
indVal_postrust_ci %>% 
  filter(assessment == "rust_ctrl") %>% 
  mutate(color0 = case_when(lowerCI == 0 ~ "nosig", TRUE ~ "sig")) %>% 
ggplot(aes(x = treatment, y = stat, color = color0)) +
  facet_grid(rows = vars(compound), cols = vars(resistance_class)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI)) +
  labs(x = "", 
       y = "Single-group indicator statistic", 
       title = "Indicator terpenes in post-rust inoculation seedlings, rust controls") +
  scale_color_manual(name = "", values = c("red", "black")) +
  theme_bw() +
  theme_bgl +
  theme(legend.position = "none")
#' 
#' ### Rust-inoculated seedlings
#' #### Indicators in QDR seedlings
#+ qdr_post_inoc
indic_post("QDR", "rust_inoc")
#' No terpenes are indicators for one treatment group. 
#' Ocimene, a_terpineol, and abietic are associated with all treatments, as opposed to controls.
#' 
#' #### Indicators in susceptible seedlings
#+ sus_post_inoc
indic_post("susceptible", "rust_inoc")
#' No terpenes are indicators for one treatment group. 
#' Ocimene, palustric, and abietic are associated with all treatments, as opposed to controls.
#' 
#' #### Indicators in MGR seedlings
#+ mgr_post_inoc
indic_post("MGR", "rust_inoc")
#' No terpenes are indicators for one treatment group. 
#' Only ocimene is associated with symbiont treatments, as opposed to controls. 
#' 
#' #### Summary
#' - **Terpenes in the post-rust, induced seedlings also often segregate along symbiont treatments vs. controls.
#' Consistent indicators in this group include ocimene and abietic acid.**
#' - **Results among resistance classes were spottier here. Where they were identified, indicators' patterns
#' were similar among resistance classes, but in several cases, indicators were only identified for one or two 
#' resistance classes.**
#' - **Note that rust inoculation is a big hammer on terpenes. An indicator analysis performed on treatment=control 
#' and rust_ctrl vs. rust_inoc assessments revealed that 21 of 26 terpenes were indicators, most of rust_inoc 
#' (not shown). This didn't seem an interesting result given the lethality of this disease.**
#' 
#' ### Plot of indicators and confidence intervals
#' The plot below shows indicator statistics and confidence intervals on single-group comparisons. 
#' The statistic shown may not match a significant pooled-group statistic if one was found using 
#' `multipatt()`. Confidence intervals are based on boostrap replication in `strassoc()` (n=1000). 
#' Confidence intervals which overlap zero mean that the statistic is non-significant. Zero-overlapping
#' CIs are shown as red on the plot.  
#+ indVal_rustinoc_plot,echo=FALSE,fig.dim=c(9,12)
indVal_postrust_ci %>% 
  filter(assessment == "rust_inoc") %>% 
  mutate(color0 = case_when(lowerCI == 0 ~ "nosig", TRUE ~ "sig")) %>% 
ggplot(aes(x = treatment, y = stat, color = color0)) +
  facet_grid(rows = vars(compound), cols = vars(resistance_class)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI)) +
  labs(x = "", 
       y = "Single-group indicator statistic", 
       title = "Indicator terpenes in post-rust inoculation seedlings, rust treated") +
  scale_color_manual(name = "", values = c("red", "black")) +
  theme_bw() +
  theme_bgl +
  theme(legend.position = "none")

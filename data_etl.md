Extract, transform, and load data
================
Beau Larkin
2022-12-08

- <a href="#description" id="toc-description">Description</a>
- <a href="#package-and-library-installation"
  id="toc-package-and-library-installation">Package and library
  installation</a>

# Description

This script pulls data from separate files in the **Database** folder,
adds a few new variables, and checks data for consistency. This script
will be called from other analysis scripts so that the ETL process isn’t
repeated, making the analysis scripts a little shorter.

# Package and library installation

Note that messages and code are often hidden in this notebook for
brevity.

``` r
# Package and library installation
packages_needed <- c("tidyverse", "magrittr", "knitr")
packages_installed <- packages_needed %in% rownames(installed.packages())
```

``` r
if (any(! packages_installed))
  install.packages(packages_needed[! packages_installed])
```

``` r
for (i in 1:length(packages_needed)) {
  library(packages_needed[i], character.only = T)
}
```

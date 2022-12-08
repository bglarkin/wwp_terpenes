Extract, transform, and load data
================
Beau Larkin
2022-12-08

- <a href="#package-and-library-installation"
  id="toc-package-and-library-installation">Package and library
  installation</a>

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

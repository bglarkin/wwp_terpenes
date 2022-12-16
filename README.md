# wwp_terpenes
Analysis of western white pine (*Pinus monitcola*) terpenes in response to blister rust, fungal endophytes, and ectomycorrhizae.

## Source data
The [Database]() directory contains raw csv files which are extracted, transformed, and loaded by [data_etl.md](data_etl.R). The ETL script is called separately in each analysis script for brevity.
- Raw script: [data_etl.R](data_etl.md)

## Global ordination and permutation tests
Ordinations of class centroids for each resistance class of seedlings. Locations based on all terpene compounds. Permutation tests conducted with `adonis2()` and 1999 permutations. 
- Report format: [terpenes_globalTests.md](terpenes_globalTests.md)
- Raw script: [terpenes_globalTests.R](terpenes_globalTests.R)

## Pairwise comparisons
1. Permutation tests are used to test pairwise multivariate differences in terpene composition
within resistance classes and assessments, and between symbiont controls and treatments. The tests are 
permuted within experimental greenhouse blocks.
2. Permutations tests as in #1, but between induced and control seedlings after the rust inoculation step
(**rust_inoc** vs. **rust_ctrl**)
- Report format: [terpenes_pairwiseTests.md](terpenes_pairwiseTests.md)
- Raw script: [terpenes_pairwiseTests.R](terpenes_pairwiseTests.R)

## Indicator terpenes 

- Report format: [terpenes_indicators.md](terpenes_indicators.md)
- Raw script: [terpenes_indicators.R](terpenes_indicators.R)

## Other files
- [gg_style.txt](gg_style.txt) contains adjustments to ggplot themes to make figures consistent within this repository.

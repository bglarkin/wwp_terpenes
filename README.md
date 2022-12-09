# wwp_terpenes
Analysis of western white pine (*Pinus monitcola*) terpenes in response to blister rust, fungal endophytes, and ectomycorrhizae.

## Source data
The [Database]() directory contains raw csv files which are extracted, transformed, and loaded by [data_etl.md](data_etl.R). The ETL script is called separately in each analysis script for brevity.
- Raw script: [data_etl.R[(data_etl.md)

## Global ordination and permutation tests
Ordinations of class centroids for each resistance class of seedlings. Locations based on all terpene compounds. Permutation tests conducted with `adonis2()` and 1999 permutations. 
- Report format: [terpenes_globalTests.md](terpenes_globalTests.md)
- Raw script: [terpenes_globalTests.R](terpenes_globalTests.R)

## Other files
- [gg_style.txt](gg_style.txt) contains adjustments to ggplot themes to make figures consistent within this repository.

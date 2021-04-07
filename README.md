# nadiv
[![](https://www.r-pkg.org/badges/version/nadiv)](https://cran.r-project.org/package=nadiv)
[![](https://cranlogs.r-pkg.org/badges/grand-total/nadiv)](https://cranlogs.r-pkg.org/badges/grand-total/nadiv)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4667837.svg)](https://doi.org/10.5281/zenodo.4667837)


R package that constructs (non)additive genetic relationship matrices, and their inverses, from a pedigree to be used in linear mixed effect models (A.K.A. the 'animal model'). Also includes other functions to facilitate the use of animal models. Some functions have been created to be used in conjunction with the R package for ASReml software.

## See the latest developments:
  - nadiv [NEWS page](https://github.com/matthewwolak/nadiv/blob/master/NEWS.md)

## Overview of main branches:
  - `master` branch is the most recent production version (typically the same as what is available from the [R CRAN mirrors](https://cran.r-project.org/))
 
  - `devel` branch is a preview of the next release which _should_ be functional and error/bug free, but proceed with caution


## To obtain nadiv:
  - From [R](https://CRAN.R-project.org/):
    - see the package page for the latest release of [nadiv on CRAN](https://CRAN.R-project.org/package=nadiv) where you can download the source.
    - install the latest release of the package directly in R:
```R
   install.packages("nadiv")
```
    - then select your favorite [CRAN mirror](https://CRAN.R-project.org/)
   
  - From GitHub:
    - clone or download the latest development version here
    - install the latest development version directly in R using the `devtools` package [https://github.com/r-lib/devtools](https://github.com/r-lib/devtools):
```R
   library(devtools); install_github("matthewwolak/nadiv")
```


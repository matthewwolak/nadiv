# nadiv
R package that constructs (non)additive genetic relationship matrices, and their inverses, from a pedigree to be used in linear mixed effect models (A.K.A. the 'animal model'). Also includes other functions to facilitate the use of animal models. Some functions have been created to be used in conjunction with the R package for ASReml software.

## To obtain nadiv:
 * From [R](http://cran.r-project.org/):
   * see the package page for the latest release of [nadiv on CRAN](http://cran.r-project.org/web/packages/nadiv/index.html) where you can download the source.
   * install the latest release of the package directly in R:
   ```R
   install.packages("nadiv")
   ```
   then selecting your favorite [CRAN mirror](http://cran.r-project.org/)
   
 * From GitHub:
   * clone or download the latest development version here
   * install the latest development version directly in R using the `devtools` package [https://github.com/hadley/devtools](https://github.com/hadley/devtools):
   ```R
   library(devtools); install_github("matthewwolak/nadiv")
   ```


## Test environments
* Ubuntu 14.04, R 3.2  #FIXME:
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'R6'

  R6 is a build-time dependency.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of nadiv 
(https://github.com/wch/checkresults/blob/master/httr/r-release). 
All packages that I could install passed except:

* 

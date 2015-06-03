## Test environments
* Ubuntu 14.04
  * R 3.2.0 (2015-04-16)
* win-builder (devel and release): http://win-builder.r-project.org/
  * R Under development (unstable) (2015-06-02 r68457)
  * R version 3.2.0 (2015-04-16)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Matthew Wolak <matthewwolak@gmail.com>'
  Possibly mis-spelled words in DESCRIPTION:
    ASReml (5:328)
    VSN (5:386)
    Suggests or Enhances not in mainstream repositories:
      asreml
  
  * spellings are, in fact, correct.
  package 'asreml' availability is noted in DESCRIPTION


* checking package dependencies ... NOTE
  Package which this enhances but not available for checking: 'asreml'

  * see response to above NOTE


## Downstream dependencies
I have also run R CMD check on downstream dependencies of nadiv: 
  dmm 
All packages that I could install passed 

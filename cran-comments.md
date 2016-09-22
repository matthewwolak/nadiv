# Test environments
* Ubuntu 14.04
  * R 3.3.1 (2016-06-21) x86_64-pc-linux-gnu (64-bit)
* win-builder (devel and release): http://win-builder.r-project.org/
  * R  Under development (unstable) (2016-09-22 r71333), platform: x86_64-w64-mingw32 (64-bit)
  * R version 3.3.1 (2016-06-21), platform: x86_64-w64-mingw32 (64-bit) 


# R CMD check results
There were no ERRORs or WARNINGs.

* There were 2 NOTEs:
  * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Matthew Wolak <matthewwolak@gmail.com>'
    Suggests or Enhances not in mainstream repositories:
    asreml

  * checking package dependencies ... NOTE
    Package which this enhances but not available for checking: 'asreml'

* No public repository is available for package 'asreml', however, availability (with web address) is noted in DESCRIPTION.


# Downstream dependencies
I have also run R CMD check on the downstream dependency of nadiv: 
  dmm, optiSel
All packages installed and passed 

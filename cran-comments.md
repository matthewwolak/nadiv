# Test environments
  - Ubuntu 17.10
    - R 3.4.3 (2017-11-30) x86_64-pc-linux-gnu (64-bit)

  - win-builder (devel and release): http://win-builder.r-project.org/
    - R version 3.5.0 (2018-04-23), platform: x86_64-w64-mingw32 (64-bit) 
    - R version 3.4.4 (2018-03-15), platform: x86_64-w64-mingw32 (64-bit) 


# R CMD check results
There were no ERRORs or WARNINGs.

  - There were 2 NOTEs:

    - checking CRAN incoming feasibility ... NOTE
        - Maintainer: 'Matthew Wolak <matthewwolak@gmail.com>'
        - Suggests or Enhances not in mainstream repositories: asreml
    - checking package dependencies ... NOTE
        - Package which this enhances but not available for checking: 'asreml'

  - No public repository is available for package 'asreml', however, availability (with web address) is noted in DESCRIPTION.

  - All other WARNINGs and NOTEs currently listed under CRAN Package Check Results have been addressed with this update. 


# Downstream dependencies
I have also run R CMD check on the downstream dependency of nadiv: 
  dmm
All packages installed and passed 

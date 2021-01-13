# Test environments
  - Ubuntu 20.04.1 LTS
    - R version 4.0.3 (2020-10-10) x86_64-pc-linux-gnu (64-bit)

  - win-builder (devel, release, and old version): https://win-builder.r-project.org/
    - R version 3.5.3 (2019-03-11), platform: x86_64-w64-mingw32 (64-bit) 
    - R version 3.6.1 (2019-07-05), platform: x86_64-w64-mingw32 (64-bit) 



# R CMD check results
There were no ERRORs or WARNINGs.

  - There were 2 NOTEs:

    - checking CRAN incoming feasibility ... NOTE
        - Maintainer: 'Matthew Wolak <matthewwolak@gmail.com>'
        - Suggests or Enhances not in mainstream repositories: asreml
    - checking package dependencies ... NOTE
        - Package which this enhances but not available for checking: 'asreml'

  - No public repository is available for package 'asreml', however, availability (with web address) is noted in DESCRIPTION.

  - All ERRORs, WARNINGs, and NOTEs currently listed under CRAN Package Check Results have been addressed with this update. 


# Downstream dependencies
I have also run R CMD check on the downstream dependency of nadiv: 
  dmm
All packages installed and passed 

# Test environments
  - Ubuntu 20.04.2 LTS
    - R version 4.0.4 (2020-02-15) x86_64-pc-linux-gnu (64-bit)

  - [win-builder](https://win-builder.r-project.org/)
    - R version 4.0.5 (2021-03-31), platform: x86_64-w64-mingw32 (64-bit)
    - R version 3.6.3 (2020-02-29), platform: x86_64-w64-mingw32 (64-bit)
    - R-devel (2021-04-05 r80144), platform: x86_64-w64-mingw32 (64-bit) 
    
  - R-hub (`devtools::check_rhub(".", interactive = FALSE)`)
    - Ubuntu Linux 20.04.1 LTS, R-release 4.0.5 (2021-03-31), GCC
    - Windows Server 2008 R2 SP1, R-devel (2021-04-01 r80137), 32/64 bit
    - Fedora Linux, R-devel (2021-04-05 r80145)
    - Debian Linux, R-devel (2020-07-31 r78945), GCC ASAN/UBSAN

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
I have also run R CMD check on the reverse dependencies and imports of nadiv: 
  dmm, optiSel
  
All packages installed and passed checks

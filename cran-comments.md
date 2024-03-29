# Test environments
  - Ubuntu 20.04.5 LTS
    - R version 4.2.2 (2022-11-10 r83330) x86_64-pc-linux-gnu (64-bit)

  - [win-builder](https://win-builder.r-project.org/)
    - R version 4.2.2 (2022-10-31), platform: x86_64-w64-mingw32 (64-bit)
    - R-devel (2022-12-06 r83409), platform: x86_64-w64-mingw32 (64-bit) 
    - R old version 4.1.3 (2022-03-10)
    
  - R-hub (`devtools::check_rhub(".", interactive = FALSE)`)
    - Ubuntu Linux 20.04.1 LTS, R-release 4.2.2 (2022-11-10), GCC
    - Windows Server 2022 R-devel (2022-10-11 r83083), 32/64 bit
    - Fedora Linux, R-devel (2022-12-06 r83409), clang, gfortran
    - Debian Linux, R-devel (2022-12-05 r83406), GCC ASAN/UBSAN

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

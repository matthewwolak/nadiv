# Test environments
  - Ubuntu 20.04.6 LTS
    - R version 4.3.2 (2023-10-31 r85441) x86_64-pc-linux-gnu

  - [win-builder](https://win-builder.r-project.org/)
    - R version 4.3.2 (2023-10-31 ucrt)
    - R Under development (unstable) (2024-01-16 r85812 ucrt)
    
  - R-hub (`devtools::check_rhub(".", interactive = FALSE)`)
    - Ubuntu 20.04.6 LTS, R version 4.3.2 (2023-10-31)
    - Windows Server 2022 x64 (build 20348), R Under development (unstable) (2023-11-18 r85554 ucrt)
    - Fedora Linux 36 (Container Image), R Under development (unstable) (2023-12-26 r85738)

  - Mac OS (`devtools::check_mac_release(".")`)
    - macOS 13.3.1 (22E261) R version 4.3.0 Patched (2023-05-18 r84451)
    
# R CMD check results
There were no ERRORs or WARNINGs.

  - There were 2 NOTEs when checking `nadiv`:

    - checking CRAN incoming feasibility ... [18s] NOTE
        - Maintainer: 'Matthew Wolak <matthewwolak@gmail.com>'
        - New submission
        - Package was archived on CRAN
        - CRAN repository db overrides: X-CRAN-Comment: Archived on 2023-12-06 as issues were not corrected in time.
        - Suggests or Enhances not in mainstream repositories:
  	    asreml
    - checking package dependencies ... NOTE
        - Package which this enhances but not available for checking: 'asreml'

  - No public repository is available for package 'asreml', however, availability (with web address) is noted in DESCRIPTION.

  - All ERRORs, WARNINGs, and NOTEs from the archived version under CRAN Package Check Results have been addressed with this update. 


# Downstream dependencies
I have also run R CMD check (`devtools::check(..., manual = FALSE, cran = TRUE, remote = TRUE, vignettes = FALSE)` on the reverse dependencies and imports of nadiv: 
  
  the latest versions archived on CRAN dmm (2.1.9), optiSel (2.0)
  
The package `dmm` installed and passed all checks.

Package `optiSel` failed on installation because the vignettes could not be built. I have been corresponding with the author of `optiSel` and will support them on any `nadiv` related issues regarding their package as they prepare a new version for CRAN.


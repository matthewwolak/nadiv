# Test environments
  - Ubuntu 20.04.6 LTS
    - R version 4.4.0 (2024-04-24 r86474) x86_64-pc-linux-gnu

  - [win-builder](https://win-builder.r-project.org/)
    - R version 4.3.3 (2024-02-29 ucrt) x86_64-w64-mingw32 (64-bit)
    - R Under development (unstable) (2024-05-16 r86559 ucrt) x86_64-w64-mingw32
    - R version 4.4.0 (2024-04-24 ucrt) x86_64-w64-mingw32
    
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

  - All ERRORs, WARNINGs, and NOTEs from the archived version under CRAN Package Check Results have been addressed with this update. All further correspondence with CRAN personnel regarding changes required after a code review have been addressed with this update.


# Downstream dependencies
I have also run R CMD check (`tools::check_packages_in_dir(..., Ncpus = 2, check_args = c("--as-cran", ""), reverse = NULL, clean = FALSE)`) on the reverse dependencies and imports of nadiv: 
  
  - the latest versions archived on CRAN that depend on `nadiv` are `dmm` (2.1-8), `optiSel` (2.0.7)
  
Both packages installed and passed all checks.



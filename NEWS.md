# 2.13.3  Released 30 May 2015
 ## NEW 
   * Added `LDtL()`, a function to take the LDL' Cholesky decomposition of a matrix (not currently exported).
   * Added `founderLine()` which traces all individuals back to either the paternal or maternal founder
   * `grfx()` now has a new argument to allow user to supply the standard normal deviates instead of generating them within the function.
     * extended the warn argument to apply to the warning when `incidence = NULL`
     * updated grfx.Rd with an example illustrating the stdnorms argument
     * added `...` argument to `drfx()` so that arguments for the internal use of `grfx()` can be supplied to `drfx()`.
   * argument now allows specified prefix for all identities in a pedigree generated from `simPedHS()` or `simPedDFC()`.
   * argument added that specifies output format of `ggcontrib()`, default is "matrix"

 ## SMALL CHANGES
   * removed 'asreml' from suggests in the package DESCRIPTION file.
   * changed `pcc()` return `FALSE` if the object (asreml) shows the log-likelihood did not converge
   * added `silent = FALSE` agrument to `pcc()` so that the default can be changed to not show messages
     * helpful in simulations where a lot of output would be printed on screen
   * changed the signs associated with likelihood ratio test statistics, etc.
     * changed signs in `proLik()` so that profile likelihoods should be "valleys" (instead of "hills", as they were in versions previous to 2.13.3)



# 2.13.2  Released 20 June 2014
 ## SMALL CHANGES
   * added the calculation of the log determinant of the A matrix to `makeAinv()`
     * log(det(A^-1)) = log(1) - log(det(A))
     * uses property of determinants that det(A^-1) = 1 / det(A) = det(A)^-1.
   * changed methods underlying `makeA()` to use inverse of cholesky factorization of the A-inverse matrix
     * `base::chol2inv()` to obtain A
       * informally seems faster unless A is dense.



# 2.13  Released 16 June 2014
 ## NEW
   * added the `prepPed()` to prepare pedigrees for use in other functions
   * exported `makeAinv()`
   * added `ggcontrib()` so that genetic group contributions can be calculated
     * still need to implement functionality for "fuzzy classification"
   * support for selfing
   * added the `pin()` and `pcc()` functions for the delta method and parameter value convergence checking, respectively, for asreml type REML models.  Also added the pin.Rd and pcc.Rd help/documentation files.
     * NOTE, `pin()` is not exported in this version (need `nadiv:::pin()` to use it)
   * added `makeDufam()`, but did not export it.
     * experimental version of makeD that first sorts individuals according to generation and then dam, and then sire.
     * sticks individuals with the same parents next to each other in the pedigree
     * haven't implemented a parallel version of this c++ code yet (or checked function for timing/memory benefits or accuracy).
   
 ## SMALL CHANGES
   * changed `makeDomEpi()` argument "Dinverse" to "invertD" to be similar to makeD()
   * added `pcc()` checks to `constrainFun()` so that only likelihood ratio test statistics of the constrained model returned if both the loglikelihood & parameter estimates have converged

# PREVIOUS (Incomplete below)
# 2.9
 ## NEW
   * enabled parallel processing (forking, so no Windows compatability) in:
     * `findDFC()`, `makeD()`, `makeDsim()`

# 2.8
   * changed name of `FindDFC()` to `findDFC()`

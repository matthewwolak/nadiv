# 2.15.0
## Small changes
   * use lower_bound algorithm for matrix lookup within c++ code
     * based on c++ <algorithm>std::lower_bound 
       * affect `makeAinv()` and `makeD()`
     * greater speedup as **A^-1** and **D** become more dense
   * create default and class 'numPed' methods for `genAssign()` and `prunePed()`
     * can greatly trim down `genAssign.numPed()` code (and to some extent `prunePed.numped()`)
     * this speeds up/uses less memory
     * since `genAssign()` and `prunePed()` are frequently called in many nadiv functions which operate on class 'numPed', this will have modest, but significant performance increases
     * thanks to [`profvis`](https://github.com/rstudio/profvis) for bringing my attention to this!

# 2.14.3 Released 20 April 2016
## NEW
  * Fuzzy classification of genetic groups to construct **A^-1**.
    * Allows individuals' phantom parents to be assigned to genetic groups with a probability. Meaning, they can be assigned to more than one genetic group.
    * To implement, the pedigree must have phantom parent identities as unique rows and a matrix of probabilities of group membership for every phantom parent in every genetic group has to be supplied to the `fuzz` argument.
    * Examples can be seen in the `makeAinv.Rd` help file or by running the following commands in `R`:

```R
?makeAinv              # launches the help documentation
example(makeAinv)      # runs the examples in the help documentation
```

    * Notably, fuzzy classification can be set to 'null', where each phantom parent is assigned to one genetic group with probability=1. This produces the same **Astar** matrix as regular genetic group methods (without fuzzy classification). See this demonstrated in the examples of the help documentation.

  * Add the `makeAstarMult()` function to create the inverse numerator relationship matrix with genetic groups (and possibly also fuzzy classification of genetic groups) through matrix multiplication instead of using direct algorithms to set this up.
    * Uses `ggcontrib()` and `makeAinv()` to create **Q** and **A^-1** directly, then multiplies these in such a way as to obtain **Astar**.
    * Examples using the two different types of pedigree formats and either with or without fuzzy classication can be seen in the `makeAstarMult.Rd` help file or run them in `R` with the command:
```R
?makeAstarMult		# launches the help documentation
example(makeAstarMult)	# runs the examples in the help documentation
```

  * Add the `F2009` dataset
    * This dataset can be used as an example for fuzzy classification of genetic groups when constructing a numerator relationship matrix with groups (i.e., with `makeAinv()`)
    * See a description in `F2009.Rd` or in R type:
```R
?F2009
```

  * Add the `simGG()` function to simulate pedigree and phenotype when immigration occurs in a focal population
    * Allows fairly fine control over a simulation. For example, the function is flexible in the: population size, number of immigrants per generation, number of generations, and both spatial and temporal trends in both focal and immigrant populations.
    * This is the function used to simulate the new `ggTutorial` dataset (below)

  * Added the `ggTutorial` dataset
    * This is a simulated dataset to be used in analyses with genetic group animal model  methods.
    * See a description in `ggTutorial.Rd` or in R type:
```R
?ggTutorial
```

  * `LRTest()` is now an exported function to do log-likelihood ratio tests

## Small changes
   * new S3 generic and methods for `makeAinv()`.
     * method dispatch is based on class of the `fuzz` argument
       * if `fuzz == NULL` then dispatch the method `makeAinv.default()`
       * if `fuzz == "matrix" | fuzz == "Matrix"` then dispatch `makeAinv.fuzzy()`

   * fix issue with `proLik()` and the confidence interval estimation
     * use `LRTest()` as basis of `constrainFun()` within `proLik()` so consistently define log-likelihood ratio test statistics
     * close issue #4 with commit [978ad610198398848d97e90c4eb57f4834a4c278](https://github.com/matthewwolak/nadiv/commits/978ad610198398848d97e90c4eb57f4834a4c278)
  

# 2.14.2 Released 5 Feb 2016
## New
   * `ggcontrib()` can now incorporate fuzzy classification of genetic groups
     * To facilitate this, the examples for `ggcontrib()` have been changed. For more information and examples, read the help documentation `ggcontrib.Rd` or in `R` type:
```R
?ggcontrib		# launches the help documentation
example(ggcontrib)	# runs the examples in the help documentation
```

## Small changes
   * fixed ordering of *f* coefficients returned by `makeAinv()`

# 2.14.1 Released 22 July 2015

# 2.14.0 Released 3 July 2015
## New
   * `makeAinv()` now can construct the augmented A-inverse matrix for genetic groups
     * This change has introduced new arguments to `makeAinv()`, however, the defaults are set to produce the normal A-inverse. For more information and examples, read the help documentation `makeAinv.Rd` or in R type:
```R
?makeAinv              # launches the help documentation
example("makeAinv")    # runs the examples in the help documentation
```
   * Improved algorithm underlying `makeAinv()` - significant speed-up!
   * Created new class `numPed` for pedigrees constructed by `numPed()`.
     * Methods for checking (`is.numPed()`) and re-ordering rows (`ronPed()`) currently available
     * To re-order the rows of an integer pedigree of class `numPed`, use `ronPed()` instead of typical subsetting operators (e.g., `'['`) to retain the class attribute `numPed`. For example:
```R
nPed <- numPed(Mrode2)
is.numPed(nPed) # TRUE
# re-order using typical R functions
nPed_sub <- nPed[order(nPed[, 2], nPed[, 3]), ]
  is.numPed(nPed_sub) # FALSE; see help via ?'[' about dropping attributes
  class(nPed_sub) # matrix
nPed_subnadiv <- ronPed(nPed, order(nPed[, 2], nPed[, 3]))
  is.numPed(nPed_subnadiv) # TRUE
  class(nPed_subnadiv) # numPed
```
   * Re-made (i.e., re-simulated) the `warcolak` dataset.
     * Codes specifying the `sex` are now `"M"` and `"F"` instead of `0` & `1`.
     * New columns added to the dataset that contain all random effects underlying the phenotype.
     * Entire code used to simulate the dataset is now an example in `warcolak.Rd`. 
   * Added two new datasets/example pedigrees: (1) `Q1988` from Quaas 1988 and (2) `Mrode3` from Mrode (2005) chapter 3. See their descriptions in `Q1988.Rd` and `Mrode3.Rd` or in R type:
```R
?Q1988
?Mrode3
```
## Small changes
   * Removed the default from the `heterogametic` argument in `makeS()`.
     * Now this has to be specified by the user!
     * For example:
```R
makeS(FG90)
``` 
won't work! Must now type:
```R
makeS(FG90, heterogametic = "0")
```
as a minimum


# 2.13.3  Released 4 June 2015

## New
   * Added `TDtT()`, a function to take the **TDT'** Cholesky decomposition of a matrix (not currently exported).
   * Added `founderLine()` which traces all individuals back to either the paternal or maternal founder
   * `grfx()` now has a new argument to allow user to supply the standard normal deviates instead of generating them within the function.
     * extended the warn argument to apply to the warning when `incidence = NULL`
     * updated grfx.Rd with an example illustrating the stdnorms argument
     * added `...` argument to `drfx()` so that arguments for the internal use of `grfx()` can be supplied to `drfx()`.
   * argument now allows specified prefix for all identities in a pedigree generated from `simPedHS()` or `simPedDFC()`.
   * argument added that specifies output format of `ggcontrib()`, default is "matrix"

## Small changes
   * removed 'asreml' from suggests in the package DESCRIPTION file.
   * changed `pcc()` return `FALSE` if the object (asreml) shows the log-likelihood did not converge
   * added `silent = FALSE` agrument to `pcc()` so that the default can be changed to not show messages
     * helpful in simulations where a lot of output would be printed on screen
   * changed the signs associated with likelihood ratio test statistics, etc.
     * changed signs in `proLik()` so that profile likelihoods should be "valleys" (instead of "hills", as they were in versions previous to 2.13.3)



# 2.13.2  Released 20 June 2014
## Small changes
   * added the calculation of the log determinant of the A matrix to `makeAinv()`
     * log(det(A^-1)) = log(1) - log(det(A))
     * uses property of determinants that det(A^-1) = 1 / det(A) = det(A)^-1.
   * changed methods underlying `makeA()` to use inverse of cholesky factorization of the A-inverse matrix
     * `base::chol2inv()` to obtain A
       * informally seems faster unless A is dense.



# 2.13  Released 16 June 2014
## New
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

## Small changes
   * changed `makeDomEpi()` argument "Dinverse" to "invertD" to be similar to makeD()
   * added `pcc()` checks to `constrainFun()` so that only likelihood ratio test statistics of the constrained model returned if both the loglikelihood & parameter estimates have converged

# PREVIOUS (Incomplete below)
# 2.9
## New
   * enabled parallel processing (forking, so no Windows compatability) in:
     * `findDFC()`, `makeD()`, `makeDsim()`

# 2.8
   * changed name of `FindDFC()` to `findDFC()`

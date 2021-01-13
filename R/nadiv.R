#' (Non)Additive Genetic Relatedness Matrices in Animal Model Analyses
#' 
#' Constructs (non)additive genetic relationship matrices, and their inverses,
#' from a pedigree to be used in linear mixed effect models (A.K.A. the 'animal
#' model'). Also includes other functions to facilitate the use of animal
#' models. Some functions have been created to be used in conjunction with the
#' R package for ASReml software, which can be obtained upon purchase from
#' VSN international (<http://www.vsni.co.uk/software/asreml>).
#' 
#' @aliases nadiv-package nadiv
#' @useDynLib nadiv, .registration = TRUE
#' @importFrom methods as is new
#' @importFrom graphics abline plot
#' @importFrom stats as.formula deriv na.omit optimize
#' @importFrom stats pchisq qchisq qnorm rnorm sd
#' @import Matrix
#' @author Matthew Wolak \email{matthewwolak@@gmail.com}
#' @examples
#' \dontrun{
#'   
#' }
"_PACKAGE"





# nadiv Cleanup: Unload DLL when library unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("nadiv", libpath)
}


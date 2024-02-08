#' Transforms ASReml-R gamma sampling variances to component scale
#' 
#' The inverse of the Average Information matrix in an ASReml-R object produces
#' the sampling variances of the (co)variance components on the gamma scale.
#' This function scales these variances to the original component scale.  This
#' allows for Confidence Intervals to be constructed about the variance
#' component estimates.
#' 
#' 
#' @param asr.object Object from a call to \code{asreml}
#' @return Returns a numeric vector of variances for each variance component in
#'   an ASReml-R model.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' 
#'   \dontrun{
#'     library(asreml)
#'     ginvA <- ainverse(warcolak)
#'     ginvD <- makeD(warcolak[, 1:3])$listDinv
#'       attr(ginvD, "rowNames") <- as.character(warcolak[, 1])
#'       attr(ginvD, "INVERSE") <- TRUE
#'     warcolak$IDD <- warcolak$ID
#'     warcolak.mod <- asreml(trait1 ~ sex,
#'      random = ~ vm(ID, ginvA) + vm(IDD, ginvD), 
#' 	data = warcolak) 
#'     summary(warcolak.mod)$varcomp
#'     sqrt(varTrans(warcolak.mod))  # sqrt() so can compare with standard errors from summary
#'    }
#' 
#' @export
varTrans <- function(asr.object){
    if(asr.object$sigma2 == 1){
       vars <- diag(aiFun(asr.object))
    } else{
       Rcomp <- which(asr.object$gammas == 1.00)
       AI <- aiFun(asr.object)
       comps <- asr.object$gammas * asr.object$sigma2
       vars <- c(((asr.object$gammas[-Rcomp]^2) * diag(AI)[Rcomp] + comps[Rcomp]^2 * diag(AI)[-Rcomp] + 2*asr.object$gammas[-Rcomp]*comps[Rcomp]*AI[Rcomp, -Rcomp]), diag(AI)[Rcomp])
      }
vars
}


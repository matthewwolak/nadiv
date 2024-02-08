#' Confidence Intervals for Variance Components
#' 
#' Produces the 1-alpha Upper and Lower Confidence Limits for the variance
#' components in an ASReml-R model.
#' 
#' Variances from the inverse of the Average Information matrix of an ASReml
#' model are translated according to the \code{\link{varTrans}} function and
#' used in constructing the 1-alpha Confidence Interval.
#' 
#' @param asr.model Object from a call to \code{asreml}
#' @param Dimnames A vector of characters if names are desired for the output.
#'   If not specified, the default labels from the \code{asreml} object will be
#'   used.
#' @param alpha A numeric value indicating the level of Type I error for
#'   constructing the Confidence Intervals.
#'
#' @return A \code{matrix} is returned with a row for each variance component. 
#'   The three columns correspond to the Lower Confidence Limit, estimate from
#'   the \code{asreml} model, and Upper Confidence Limit for each variance 
#'   component.
#' @note The vector of \code{Dimnames} should match the same order of variance
#' components specified in the model.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{aiFun}} \code{\link{proLik}}
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
#'     aiCI(warcolak.mod)    
#'    }
#' 
#' @export
aiCI <- function(asr.model, Dimnames = NULL, alpha = 0.05)
{
   za2 <- qnorm(alpha/2, mean = 0, sd = 1) 
   hii.vec <- varTrans(asr.model)
   theta.vec <- asr.model$gammas * asr.model$sigma2
   UCL <- theta.vec - za2*sqrt(hii.vec)
   LCL <- theta.vec + za2*sqrt(hii.vec)
  CIframe <- cbind(LCL, theta.vec, UCL)
  if(!is.null(Dimnames)) dimnames(CIframe)[[1]] <- Dimnames
  dimnames(CIframe)[[2]] <- c("LCL", "estimate", "UCL")
return(CIframe)
}


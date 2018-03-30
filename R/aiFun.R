#' Sampling (co)variances
#' 
#' This function returns the sampling (co)variances of the variance components
#' fitted in an mixed model solved using the Average Information algorithm
#' 
#' The inverse of the Average Information matrix provides the sampling
#' (co)variance of each (co)variance component in the random portion of the
#' mixed model. If a model from the ASReml-R function is supplied (\code{model}
#' is not NULL), this function extracts the inverse of the AI matrix from an
#' ASReml-R model and organizes it so that the sampling covariances between
#' random terms are the off-diagonals and the sampling variances of random
#' terms are located along the diagonal. The order of the variances along the
#' diagonal is the same as the order entered in the random section of the
#' \code{asreml} function. This is also the same order as the rows of a call to
#' the summary function, \code{summary(model)$varcomp}.
#' 
#' If \code{model} is NULL then \code{AI.vec} should contain the vector of
#' values from an Average Information matrix. The function will then
#' reconstruct this matrix, invert it, and supply the sampling (co) variances
#' for the random terms in the model as described above. Note, either
#' \code{model} or \code{AI.vec} must be supplied, but not both.
#' 
#' @usage aiFun(model = NULL, AI.vec = NULL, inverse = TRUE, Dimnames = NULL)
#' @param model A model object returned by a call to the \code{asreml}
#' function.
#' @param AI.vec A numeric vector of the Average Information matrix. The order
#' must be the row-wise lower triangle of the matrix (including the diagonal).
#' @param inverse A logical indicating whether the elements of the
#' \emph{inverse} Average Information matrix are being provided. If FALSE, the
#' Average Information matrix (and not its inverse) is being supplied.
#' @param Dimnames A vector of characters if names are desired for the output
#' (co)variance matrix. If not specified, either the default labels from the
#' \code{asreml} object will be used or the rows and columns will be
#' un-labeled.
#' @return A matrix of k x k dimensions is returned, if k is the number of
#' (co)variance components estimated in the model. Sampling covariances are
#' above and below the diagonal while variances are located along the diagonal.
#' If \code{Dimnames} is specified, the row and column names are assigned
#' according the vector of names in this argument.
#' @note The vector of \code{Dimnames} should match the same order of variance
#' components specified in the model.
#' @author \email{matthewwolak@@gmail.com}
#' @references Gilmour, A.R., Gogel, B.J., Cullis, B.R., & Thompson, R. 2009.
#' ASReml User Guide Release 3.0. VSN International Ltd., Hemel Hempstead, UK.
#' @examples
#' 
#'   \dontrun{
#'     library(asreml)
#'     ginvA <- asreml.Ainverse(warcolak)$ginv
#'     ginvD <- makeD(warcolak[,1:3])$listDinv
#'     warcolak$IDD <- warcolak$ID
#'     warcolak.mod <- asreml(trait1 ~ sex, random = ~ped(ID) + giv(IDD), 
#' 	ginverse = list(ID = ginvA, IDD = ginvD), data = warcolak) 
#'     summary(warcolak.mod)$varcomp
#'     aiFun(model = warcolak.mod, Dimnames = c("Va", "Vd", "Ve"), inverse = TRUE)    
#'    }
#' 
#'   output <- c(7.3075921, 7.0635161, 12.3423380, 1.9539486, 2.7586340, 0.6626111)
#'   aiFun(AI.vec = output, inverse = FALSE, Dimnames = c("Va", "Vd", "Ve"))
#' 
#' @export aiFun
aiFun <- function(model = NULL, AI.vec = NULL, inverse = TRUE, Dimnames=NULL)
{
  if(!is.null(model)){
    AI.vec <- model$ai
  }

  dimAI <- sqrt(length(AI.vec) * 2 + 0.25) - 0.5
  AI <- matrix(0, dimAI, dimAI)
  AI[which(upper.tri(AI, diag = TRUE) == TRUE)] <- AI.vec
  AI[which(lower.tri(AI) == TRUE)]<-t(AI)[which(lower.tri(AI) == TRUE)]
  if(inverse == FALSE) AI <- solve(AI)    

  if(is.null(Dimnames)){Dimnames <- names(model$gammas)}
  dimnames(AI) <- list(Dimnames, Dimnames)
       
AI
}


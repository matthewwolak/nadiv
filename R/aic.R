#' Akaike Information Criterion
#' 
#' Calculates AIC/AICc values, AIC differences, Likelihood of models, and model
#' probabilities.
#' 
#' Calculations and notation follows chapter 2 of Burnham and Anderson (2002).
#' 
#' @param logLik A vector of model log-Likelihoods
#' @param fp A vector containing the numbers of free parameters of each model
#'   included in the logLik vector
#' @param n An optional vector of sample sizes for each model.  Used to
#'   calculate AICc (small sample un-biased AIC).
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{AIC }{vector containing AIC/AICc (depending on value of \code{n})}
#'     \item{delta_AIC }{vector containing AIC differences from the minimum 
#'       AIC(c)}
#'     \item{AIClik }{vector containing likelihoods for each model, given the 
#'       data.  Represents the relative strength of evidence for each model.}
#'     \item{w }{Akaike weights.}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @references Burnham, K.P. and D.R. Anderson. 2002. Model Selection and
#'   Multimodel Inference. A Practical Information-Theoretic Approach, 2nd edn.
#'   Springer, New York.
#' @examples
#' 
#'    aic(c(-3139.076, -3136.784, -3140.879, -3152.432), c(8, 7, 8, 5)) 
#' 
#' @export
aic <- function(logLik, fp, n = NULL){
   if(is.numeric(n)){
      AICc <- -2*(logLik - fp * (n / (n - fp - 1)))
      delta_AIC <- AICc - min(AICc)
   } else{
        AIC <- -2*(logLik - fp)
        delta_AIC <- AIC - min(AIC)
     } 
   AIClik <- exp(-0.5*delta_AIC)
   w <- AIClik / sum(AIClik)
   
   if(is.numeric(n)){
      return(list(AICc = AICc, delta_AIC = delta_AIC, AIClik = AIClik, w = w))
   } else{
        return(list(AIC = AIC, delta_AIC = delta_AIC, AIClik = AIClik, w = w))
     }
}


# Likelihood Ratio Test:
# Tests the hypothesis that the reduced model offers a better fit
# Helpful reading: section 6.4.1.1 of Bolker 2008. pp. 189-194
#NOTE: sometimes ASReml-R returns positive log-likelihoods, other times negative log-likelihoods. My experience so far is that positive values are returned when there are often boundary parameters or the model is having trouble fitting the data. 
#	 ***!!!ASSUMPTION that could be implemented!!!*** 
#	when both are positive, the likelihood is being maximized (i.e., the greater value is a greater likelihood.
LRTest <- function(full, reduced, df = 1, boundaryCorrection = FALSE){
    if(sign(full) != sign(reduced)){
       stop("Signs of the log-likelihoods are opposite - or 1 log-likelihood is zero...don't know what to do")
    }
    if(sign(full) > 0) stop("positive log-likelihoods")
    # negative log-likelihoods: better fit has greater log-likelihood (less negative)
    if(sign(full) < 0){
       lambda <- 2*(full - reduced)
    }

    if(boundaryCorrection & df == 1){
       lrtP <- 0.5*(pchisq(lambda, df = df, lower.tail = FALSE))
    } else{
         lrtP <- pchisq(lambda, df = df, lower.tail = FALSE)
      }
         
 return(list(lambda = lambda, Pval = lrtP, corrected.Pval = boundaryCorrection))
}


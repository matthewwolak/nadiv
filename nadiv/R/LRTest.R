# Likelihood Ratio Test:
# Tests the hypothesis that the reduced model offers a better fit
# Helpful reading: section 6.4.1.1 of Bolker 2008. pp. 189-194
#NOTE: sometimes ASReml-R returns positive log-likelihoods, other times negative log-likelihoods. My experience so far is that positive values are returned when there are often boundary parameters or the model is having trouble fitting the data. Regardless of how it occurs, to account for this, I will always use the absolute values of the log-likelihoods:
#	 ***!!!ASSUMPTIONs!!!*** 
#	when both are negative, the likelihood is being minimized (i.e., the more negative value is a greater likelihood)
#	when both are positive, the likelihood is being maximized (i.e., the greater value is a greater likelihood.
LRTest <- function(full, reduced, df = 1, boundaryCorrection = FALSE){
    if(sign(full) != sign(reduced)){
       stop("Signs of the log-likelihoods are opposite - or 1 log-likelihood is zero...don't know what to do")
    }
    lambda <- 2*(abs(full) - abs(reduced))
    if(boundaryCorrection & df == 1){
       lrtP <- 0.5*(pchisq(lambda, df = df, lower.tail = FALSE))
    } else{
         lrtP <- pchisq(lambda, df = df, lower.tail = FALSE)
      }
         
 return(list(lambda = lambda, Pval = lrtP, corrected.Pval = boundaryCorrection))
}


LRTest <- function(full, reduced, df = 1, boundaryCorrection = FALSE){
    lambda <- (2*(full - reduced))
    if(boundaryCorrection & df == 1){
       lrtP <- 0.5*(pchisq(lambda, df = df, lower.tail = FALSE))
    } else{
         lrtP <- pchisq(lambda, df = df, lower.tail = FALSE)
      }
         
 return(list(lambda = lambda, Pval = lrtP, corrected.Pval = boundaryCorrection))
}


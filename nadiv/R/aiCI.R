aiCI <- function(asr.model, Dimnames = NULL, alpha = 0.05)
{
   za2 <- qnorm(alpha/2, mean = 0, sd = 1) 
   hii.vec <- varTrans(asr.model)
   theta.vec <- summary(asr.model)$varcomp$component
   UCL <- theta.vec - za2*sqrt(hii.vec)
   LCL <- theta.vec + za2*sqrt(hii.vec)
  CIframe <- cbind(LCL, theta.vec, UCL)
  if(!is.null(Dimnames)) dimnames(CIframe)[[1]] <- Dimnames
  dimnames(CIframe)[[2]] <- c("LCL", "estimate", "UCL")
return(CIframe)
}


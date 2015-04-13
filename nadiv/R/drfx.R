drfx <- function(G, fac, dataf){
    dataf[, fac] <- as.factor(dataf[, fac])
    d <- nrow(G)
    if(all(G == G[1,1]) & d > 1) warning("variance-covariance matrix 'G' may have caused 'chol.default(G)' error.  If so, consider subtracting 0.0001 from the covariances to make correlations < 1 or >-1")

   Z <- sparse.model.matrix(as.formula(paste("~", fac, " - 1", sep = "")), dataf)
   M <- suppressWarnings(grfx(n = ncol(Z), G = G))
   fx <- sapply(seq.int(d), FUN = function(c){ (Z %*% M[, c])@x})
 return(list(fx = fx, Z = Z))
}


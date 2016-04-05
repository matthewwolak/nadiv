TDtT <- function(A, ...){
    ch <- t(chol(A))
    dd <- diag(ch)
 return(list(T = drop0(zapsmall(ch / dd, ...)), D = Diagonal(nrow(A), dd^2)))
}


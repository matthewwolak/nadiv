TDtT <- function(A){
    ch <- t(chol(A))
    dd <- diag(ch)
 return(list(T = ch / dd, D = Diagonal(nrow(A), dd^2)))
}


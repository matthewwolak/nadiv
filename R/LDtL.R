LDtL <- function(A){
    ch <- chol(A)
    dd <- diag(ch)
 return(list(tL = ch / dd, D = Diagonal(nrow(A), dd^2)))
}


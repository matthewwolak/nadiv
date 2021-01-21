#####################################
#adapted from code written by
#Jarrod Hadfield in the 
#MCMCglmm package
######################################


#' Converts a sparse matrix into a three column format.
#' 
#' From a sparse matrix object, the three column, row ordered lower triangle of
#' non-zero elements is created.  Mostly used within other functions (i.e.,
#' \code{makeD})
#' 
#' The sparse matrix and three column format must fit CERTAIN assumptions about
#' row/column sorting and lower/upper triangle matrix.
#' 
#' Adapted from a function in the \code{MCMCglmm} package
#' 
#' @param A a sparse matrix
#' @param rownames a list of rownames from the 'A' matrix.
#' @param colnames the columns will be labeled however they are entered in
#'   this character vector
#' @return returns the list form of the sparse matrix as a \code{data.frame}
#' @seealso \code{\link[MCMCglmm]{MCMCglmm}}
#' @export
sm2list<-function(A, rownames = NULL, colnames=c("row", "column", "A"))
{
    ginv <- data.frame(Row = rep(1:length(A@p[-1]),
        diff(A@p)), Column = A@i + 1, Ainverse = A@x)
    ginv <- ginv[order(ginv$Row), ]
    ginv <- ginv[which(ginv$Row >= ginv$Column), ]
    attr(ginv, "rowNames") <- rownames
    names(ginv)<-colnames
    return(ginv)
}


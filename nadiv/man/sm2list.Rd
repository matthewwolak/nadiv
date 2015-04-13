\name{sm2list}
\alias{sm2list}
\title{Converts a sparse matrix into a three column format.
}
\description{From a sparse matrix object, the three column, row ordered lower triangle of non-zero elements is created.  Mostly used within other functions (i.e., \code{makeD})
}
\usage{
sm2list(A, rownames = NULL, colnames = c("row", "column", "A"))
}
\arguments{
  \item{A }{
a sparse matrix
}
  \item{rownames }{
a list of rownames from the 'A' matrix.
}
  \item{colnames }{
the columns will be labelled however they are entered in this character vector}
}
\details{
The sparse matrix and three column format must fit CERTAIN assumptions about row/column sorting and lower/upper triangle matrix.

Adapted from a function in the \code{MCMCglmm} package}

\value{
returns the list form of the sparse matrix as a \code{data.frame}}

\seealso{\code{\link[MCMCglmm]{MCMCglmm}}}


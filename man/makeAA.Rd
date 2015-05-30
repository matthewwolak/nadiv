\name{makeAA}
\alias{makeAA}
\title{Creates the additive by additive epistatic genetic relationship matrix
}
\description{
Given a pedigree, the matrix of additive by additive genetic relatedness (AA) among all individuals in the pedigree is returned.}
\usage{
makeAA(pedigree)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}
}
\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

The function first estimates the A matrix using \code{\link{makeA}}, then it calculates the Hadamard (element-wise) product of the A matrix with itself (A # A).}

\value{
  \item{AA }{the AA matrix in sparse matrix form}
  \item{logDet }{the log determinant of the AA matrix}
  \item{AAinv }{the inverse of the AA matrix in sparse matrix form}
  \item{listAAinv }{the three column form of the non-zero elements for the inverse of the AA matrix}
}

\author{\email{matthewwolak@gmail.com}
}


\seealso{
\code{\link{makeA}}
}

\examples{
  makeAA(Mrode2)
}


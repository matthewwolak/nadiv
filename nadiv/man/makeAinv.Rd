\name{makeAinv}
\alias{makeAinv}
\title{Creates the Inverse additive genetic relationship matrix
}
\description{This returns the inverse of the additive relationship matrix.  It can also be used to obtain coefficients of inbreeding for the pedigreed population.
}
\usage{
makeAinv(pedigree, det = FALSE)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}
  \item{det }{
Logical, indicating if the (log) determinant of the A matrix should be returned
}
}

\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

The function implements an adaptation of the Meuwissen and Luo (1992) algorithm and the code was borrowed heavily from the \code{inverseA} function by Jarrod Hadfield in the \code{MCMCglmm} package. 
}

\value{
  \item{Ainv }{the inverse of the additive genetic relationship matrix in sparse matrix form}
  \item{listAinv }{the three column form of the non-zero elements for the inverse of the additive genetic relationship matrix}
  \item{f }{the inbreeding coefficients for each individual in the pedigree}
  \item{logDet }{the log determinant of the A matrix}
}

\references{
Meuwissen, T.H.E & Luo, Z. 1992. Computing inbreeding coefficients in large populations. Genetics, Selection, Evolution. 24:305-313.
}

\author{\email{matthewwolak@gmail.com}}

\seealso{\code{\link{makeA}}}

\examples{
 makeAinv(Mrode2)
}


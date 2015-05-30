\name{makeD}
\alias{makeD}
\title{Creates the dominance genetic realationship matrix
}
\description{Given a pedigree, the matrix of coefficients of fraternity are returned - the D matrix.  Note, no inbreeding must be assumed.  Will return the inverse of the D matrix by default, otherwise this operation can be skipped if desired.
}
\usage{
makeD(pedigree, parallel = FALSE, ncores = getOption("mc.cores", 2L),
      invertD = TRUE, returnA = FALSE, det = FALSE)
}
\arguments{
  \item{pedigree }{
A pedigree with columns organized: ID, Dam, Sire}

  \item{parallel }{
Logical, indicating whether computation should be run on multiple processors at once.  See details for considerations.}

  \item{ncores }{
Number of cores to use when parallel = TRUE.  Default is maximum available.  Otherwise, set with an integer.  See details for considerations.}

  \item{invertD }{
A logical indicating whether or not to invert the D matrix}

  \item{returnA }{
Logical, indicating if the numerator relationship matrix (A) should be stored and returned.}

  \item{det }{
Logical, indicating if the determinant of the D matrix should be returned.}

}
\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

There exists no convenient method of obtaining the inverse of the dominance genetic relatedness matrix (or the D matrix itself) directly from a pedigree (such as for the inverse of A, i.e., Quaas (1995)).  Therefore, this function computes the coefficient of fraternity (Lynch and Walsh, 1998) for every individual in the pedigree with a non-zero additive genetic relatedness.  Note, the construction of the D matrix is more computationally demanding (in computing time and space requirements) than is the construction of A.

To overcome the computational difficulties, this function can enable parallel processing (see package \code{parallel} included in the R distribution) to speed up the execution.  Note this may not be possible on Windows (See \code{parallel} documentation for further information), therefore \code{parallel} = TRUE should only be used on Linux or Mac operating systems (i.e., not Windows).  The default is to use the maximum available cpus to the machine, but this can be restricted by indicating the number desired in the argument \code{ncores}.  Setting up the multi-processing takes some overhead, so no real advantage is gained for small pedigrees.  Also, since all processes are sharing a fixed amount of RAM, for very large pedigrees using many processes in parallel may not be feasible due to RAM restrictions (i.e., if each process needs "n" amount of RAM to run, then \code{ncores} should be set to = total RAM/n).  Otherwise the machine can become overworked.  

Note, for very large pedigrees \code{returnA} should be set to FALSE to avoid drastically increasing the memory requirements while making D.  When this occurs, 'NULL' is returned for the element of 'A' in the output of \code{makeD}.
}

\value{
  \item{A }{the A matrix in sparse matrix form}
  \item{D }{the D matrix in sparse matrix form}
  \item{logDet }{the log determinant of the D matrix}
  \item{Dinv }{the inverse of the D matrix in sparse matrix form}
  \item{listDinv }{the three column form of the non-zero elements for the inverse of the D matrix}
}
\references{
Quaas, R.L. 1995. Fx algorithms. An unpublished note.

Lynch M., & Walsh, B. 1998. Genetics and Analysis of Quantitative Traits. Sinauer, Sunderland, Massachusetts.}
\author{\email{matthewwolak@gmail.com}
}

\seealso{
  \code{\link{makeDsim}} 
}


\examples{
  DinvMat <- makeD(Mrode9, parallel = FALSE)$Dinv
}


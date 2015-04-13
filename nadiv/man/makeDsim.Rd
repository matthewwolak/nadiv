\name{makeDsim}
\alias{makeDsim}
\title{
Creates the dominance genetic relationship matrix through an iterative (simulation) process}
\description{
Alleles are explicitly traced through a pedigree to obtain coefficients of fraternity between pairs of individuals (the probability of sharing both alleles identical by descent).  This is accomplished in an iterative process to account for the various routes by which an allele will progress through a pedigree due to Mendelian sampling.  This is an implementation of the simulation approach of Ovaskainen et al. (2008).}
\usage{
makeDsim(pedigree, N, parallel = FALSE, ncores = getOption("mc.cores", 2L),
         invertD = TRUE, calcSE = FALSE, returnA = FALSE)
}

\arguments{
  \item{pedigree }{
A pedigree with columns organized: ID, Dam, Sire}
  \item{N }{
The number of times to iteratively trace alleles through the pedigree}
  \item{parallel }{
A logical indicating whether or not to use parallel processing.  Note, this may only be available for Mac and Linux operating systems.}
  \item{ncores }{
The number of cpus to use when constructing the dominance relatedness matrix.  Default is all available.}
  \item{invertD }{
A logical indicating whether or not to invert the D matrix}
  \item{calcSE }{
A logical indicating whether or not the standard errors for each coefficient of fraternity should be calculated}
  \item{returnA }{
Logical, indicating if the numerator relationship matrix (A) should be stored and returned.}
}

\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

\code{parallel} = TRUE should only be used on Linux or Mac operating systems (i.e., not Windows).

Ovaskainen et al. (2008) indicated that the method of calculating the D matrix (see \code{\link{makeD}}) is only an approximation.  They proposed a simulation method that is implemented here.  This should be more appropriate, especially when inbreeding occurs in the pedigree. 

The value, \code{listDsim} will list both the approximate values (returned from \code{\link{makeD}}) as well as the simulated values.  If \code{calcSE} is TRUE, these values will be listed in \code{listDsim}.
}

\value{
  \item{A }{the A matrix in sparse matrix form}
  \item{D }{the approximate D matrix in sparse matrix form}
  \item{logDetD }{the log determinant of the approximate D matrix}
  \item{Dinv }{the inverse of the approximate D matrix in sparse matrix form}
  \item{listDinv }{the three column form of the non-zero elements for the inverse of the approximate D matrix}
  \item{Dsim }{the simulated D matrix in sparse matrix form}
  \item{logDetDsim }{the log determinant of the simulated D matrix}
  \item{Dsiminv }{the inverse of the simulated D matrix in sparse matrix form}
  \item{listDsim }{the three column form of the non-zero and non-self elements for the simulated D matrix}
  \item{listDsiminv }{the three column form of the non-zero elements for the inverse of the simulated D matrix}

}
\references{
Ovaskainen, O., Cano, J.M., & Merila, J. 2008. A Bayesian framework for comparative quantitative genetics. Proceedings of the Royal Society B 275, 669-678.
}
\author{\email{matthewwolak@gmail.com}
}

\note{
This simulation can take a long time for large values of \code{N}.  If unsure, it is advisable to start with a lower N and gradually increase to gain a sense of the time required to execute a desired \code{N}.}

\seealso{
 \code{\link{makeD}} 
}


\examples{
  simDinv <- makeDsim(Mrode9, N = 1000, parallel = FALSE, invertD = TRUE,
		calcSE = TRUE)$listDsim
}


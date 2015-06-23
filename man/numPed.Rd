\name{numPed}
  \alias{numPed}
  \alias{is.numPed}
  \alias{ronPed}

\title{Integer Format Pedigree}

\description{
Conversion, checking, and row re-ordering of a pedigree in integer form of class \sQuote{numPed}.}

\usage{
numPed(pedigree, check = TRUE)
is.numPed(x)
ronPed(x, i, \dots)
}
\arguments{
  \item{pedigree }{
A three column pedigree object, where the columns correspond to: ID, Dam, & Sire}
  \item{check }{
A logical argument indicating if checks on the validity of the pedigree structure should be made, but see Details}
  \item{x }{
A pedigree of class \sQuote{\code{numPed}}}
  \item{i, \dots }{
Index specifying elements to extract or replace: see \code{\link[base]{[}}}


}
\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', '-998', or '*'.

Individuals must appear in the ID column in rows preceeding where they appear in either the Dam or Sire column. See the \code{\link[nadiv]{prepPed}} function if this is not the case.

If pedigree inherits the class "numPed" (from a previous call to \code{numPed()}) and \code{check = TRUE}, the checks are skipped. If \code{check = FALSE} any pedigree will be transformed into a pedigree consisting of integers and missing values denoted by '-998'.

Based on code from the \code{MCMCglmm} package}

\seealso{\code{\link[nadiv]{prepPed}}, \code{\link[MCMCglmm]{MCMCglmm}}, \code{\link[base]{[}}}


\value{
An S3 object of class "numPed" representing the pedigree, where individuals are now numbered from 1 to n and unknown parents are assigned a value of '-998'.
}

\examples{
(nPed <- numPed(Mrode2))
is.numPed(nPed)

# re-order and retain class 'numPed'
ronPed(nPed, order(nPed[, 2], nPed[, 3]))

}


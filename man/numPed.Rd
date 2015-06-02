\name{numPed}
\alias{numPed}
\title{Creates an integer format for a pedigree
}
\description{
Converts a pedigree into a standard integer form used by many functions in the 'nadiv' package.}
\usage{
numPed(pedigree, check = TRUE)
}
\arguments{
  \item{pedigree }{
A three column pedigree object, where the columns correspond to: ID, Dam, & Sire}
  \item{check }{
A logical argument indicating if checks on the validity of the pedigree structure should be made}

}
\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0',, '-998', or '*'.

Individuals must appear in the ID column in rows preceeding where they appear in either the Dam or Sire column.  See the \code{\link[nadiv]{prepPed}} function if this is not the case.  If \code{check = FALSE} any pedigree will be transformed into a pedigree consisting of integers and missing values denoted by '-998'.

Based on code from the \code{MCMCglmm} package}

\seealso{\code{\link[nadiv]{prepPed}}, \code{\link[MCMCglmm]{MCMCglmm}}}


\value{
An S3 object of class "numPed" representing the pedigree, where individuals are now numbered from 1 to n and unknown parents are assigned a value of '-998'.
}



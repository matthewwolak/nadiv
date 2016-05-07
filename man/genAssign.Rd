\name{genAssign}
\alias{genAssign}
\alias{genAssign.default}
\alias{genAssign.numPed}

\title{
Generation assignment}
\description{
Given a pedigree, the function assigns the generation number to which each individual belongs.}

\usage{
genAssign(pedigree, \dots)
\method{genAssign}{default}(pedigree, \dots)
\method{genAssign}{numPed}(pedigree, \dots)
}
\arguments{
  \item{pedigree }{A pedigree where the columns are ordered ID, Dam, Sire}
  \item{\dots }{Arguments to be passed to methods}
}
\details{
0 is the base population.  

Migrants, or any individuals where both parens are unknown, are assigned to generation zero.  If parents of an individual are from two different generations (e.g., dam = 0 and sire = 1), the individual is assigned to the generation following the greater of the two parents (e.g., 2 in this example).
} 


\value{
A vector of values is returned.  This vector is in the same order as the ID column of the pedigree.
}

\author{\email{matthewwolak@gmail.com}}




\name{genAssign}
\alias{genAssign}

\title{
Generation assignment}
\description{
Given a pedigree, the function assigns the generation number to which each individual belongs.}

\usage{
genAssign(pedigree)
}
\arguments{
  \item{pedigree }{A pedigree with columns organized: ID, Dam, Sire
}
}
\details{
0 is the base population.  

Migrants, or any individuals where both parens are unknown, are assigned to generation zero.  If parents of an individual are from two different generations (e.g., mom = 0 and dad = 1), the individual is assigned to the generation following the greater of the two parents (e.g., 2 in this example).
} 


\value{
A vector of values is returned.  This vector is in the same order as the ID column of the pedigree.
}



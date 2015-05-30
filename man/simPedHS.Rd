\name{simPedHS}
\alias{simPedHS}
\title{
Half-sib pedigree construction
}
\description{
Simulates a pedigree for a half-sib mating design (sometimes also called the North Carolina Design 1).
}
\usage{
simPedHS(s, d, n, uniqueDname = TRUE, prefix = NULL)
}
\arguments{
  \item{s }{
Number of sires
}
  \item{d }{
Number of dams per sire
}
  \item{n }{
Number of offspring per mating (must be > or = 2)
}
  \item{uniqueDname }{
Logical indicating if dams should have unique names within sire families or throughout the entire pedigree
}
  \item{prefix }{
Optional prefix to every identity
}
 
}
\details{
\code{n} must be greater than or equal to 2, because one male and one female offspring are produced from each mating

Some functions/calculations get bogged down if no two dams have the same ID in the entire pedigree (e.g., \code{aov}).  However, other functions must have unique identifiers for every individual.
}

\value{
A \code{data.frame} with columns corresponding to: id, dam, sire, and sex.  Sex is "M" for males and "F" for females.
}

\author{\email{matthewwolak@gmail.com}}

\seealso{
  \code{\link{simPedDFC}}
}
\examples{
  simPedHS(s = 1, d = 3, n = 2)
}


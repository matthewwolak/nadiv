\name{simPedDFC}
\alias{simPedDFC}
\title{
Double first cousin pedigree construction
}
\description{
Simulates a pedigree for the \dQuote{double first cousin} mating design (Fairbairn and Roff 2006).
}
\usage{
simPedDFC(F, gpn = 4, fsn = 4, s = 2, prefix = NULL)
}
\arguments{
  \item{F }{
Number of blocks for the design
}
  \item{gpn }{
Number of grandparents in the first/GP generation (must be >= 2)
}
  \item{fsn }{
Number of offspring in the full-sib families of the second/P generation (must be an even number >= 4)
}
  \item{s }{
Number of sires per full-sib family in the second/P generation (must be >=2)
}
  \item{prefix }{
Optional prefix to add to every identity
}

}
\details{
This is an adaption to a half-sib breeding design which also produces first cousins and double first cousins.  Double first cousins are produced by mating two brothers to two sisters (the offspring of the resulting two families are double first cousins with one another).  This is described in Fairbairn and Roff (2006) as being particularly effective for separating autosomal additive genetic variance from sex chromosomal additive genetic variance.  It is also amenable to estimating dominance variance, however, it still has difficulty separating dominance variance from common maternal environmental variance (Meyer 2008).
}
\value{
A \code{data.frame} with columns corresponding to: id, dam, sire, and sex.  Sex is \code{M} for males and \code{F} for females.
}
\references{
Fairbairn, D.J. and D.A. Roff. 2006.  The quantitative genetics of sexual dimorphism: assessing the importance of sex-linkage. Heredity 97:319-328.
 
Meyer, K. 2008.  Likelihood calculations to evaluate experimental designs to estimate genetic variances. Heredity 101:212-221. 
}

\author{\email{matthewwolak@gmail.com}}



\seealso{
  \code{\link{simPedHS}}, \code{\link{warcolak}}
}
\examples{
  DFC1 <- simPedDFC(F = 1, gpn = 2, fsn = 4, s = 2)
}


\name{warcolak}
\alias{warcolak}
\docType{data}
\title{
Pedigree and phenotypic values for a mythical population of Warcolaks}
\description{
A two trait example pedigree from the three generation breeding design of Fairbairn & Roff (2006) with two un-correlated traits.}

\usage{data(warcolak)}
\format{
  A data frame with 5400 observations on the following 5 variables.
  \describe{
    \item{\code{ID }}{a factor specifying 5400 unique individual IDs}
    \item{\code{Dam }}{a factor specifying the unique ID for each dam}
    \item{\code{Sire }}{a factor specifying the unique ID for each sire}
    \item{\code{sex} }{a vector specifying '0' if the individual is a male and '1' if it is a female}
    \item{\code{trait1} }{a numeric vector of phenotypic values}
    \item{\code{trait2} }{a numeric vector of phenotypic values}
  }
}

\details{
Unique sets of relatives are specified (Fairbairn & Roff, 2006) for a three generation breeding design.  Each set contains 72 individuals.  This pedigree reflects an experiment which produces 75 of these basic sets from Fairbairn & Roff's design.

The dataset was simulated to have two un-correlated traits with different genetic architectures.  The additive genetic, dominance genetic, and environmental (or residual) variances for both \code{trait1} and \code{trait2} are specified to be 0.4, 0.3, & 0.3, respectively.  However, the additive genetic variance for \code{trait2} can be further decomposed to autosomal additive genetic variance (0.3) and X-linked additive genetic variance (0.1).  These additive genetic effects were drawn from multivariate random normal distributions [e.g., autosomal additive genetic effects: N ~ (0, A * Va)] with means of zero and variances equal to the product of the expected variance and the relatedness (or incidence) matrix.  Because of Monte Carlo error, the actual variance in random effects will vary slightly from the amount specified in the simulation. 
}

\references{
Fairbairn, D.J. & Roff, D.A. 2006. The quantitative genetics of sexual dimorphism: assessing the importance of sex-linkage. Heredity 97, 319-328.}


\keyword{datasets}

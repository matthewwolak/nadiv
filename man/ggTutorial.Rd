\name{ggTutorial}
\alias{ggTutorial}
\docType{data}
\title{
Simulated dataset used to analyze data with genetic group animal models
}
\description{
The dataset was simulated using the \code{\link{simGG}} function so that the pedigree contains a base population comprised of founders and non-founder immigrants. These data are then used in the main manuscript and tutorials accompanying Wolak & Reid (unpublished manuscript).
}
\usage{data("ggTutorial")}
\format{
  A data frame with 6000 observations on the following 10 variables.
  \describe{
    \item{\code{id }}{an integer vector specifying the 6000 unique individual identities}
    \item{\code{dam }}{an integer vector specifying the unique dam for each individual}
    \item{\code{sire }}{an integer vector specifying the unique sire for each individual}
    \item{\code{parAvgU }}{a numeric vector of each individual's parents' average autosomal additive genetic effects}
    \item{\code{mendel }}{a numeric vector of each individual's deviation from their \code{parAvgU} autosomal additive genetic effects}
    \item{\code{u }}{a numeric vector of the total autosomal additive genetic effects underlying \code{p}}
    \item{\code{r }}{a numeric vector of the residual (environmental) effects underlying \code{p}}
    \item{\code{p }}{a numeric vector of phenotypic values}
    \item{\code{is }}{an integer vector with \code{0} for individuals born in the focal population and \code{1} for individuals born outside of the focal population, but immigrated}
    \item{\code{gen }}{an integer vector specifying the generation in which each individual was born}
  }
}

\details{
The dataset was simulated as described in the \sQuote{examples} section using the \code{\link{simGG}} function. Full details of the function and dataset can be found in Wolak & Reid (unpublished manuscript). 
}

\source{
Wolak, M.E. & J.M. Reid. Unpublished manuscript. Accounting for genetic differences among unknown parents in microevolutionary studies: How to include genetic groups in quantitative genetic animal models.
}

\examples{
 \dontrun{
  rm(list = ls())
  set.seed(102)     #<-- seed value used originally
  library(nadiv)
  # create data using `simGG()`
  ggTutorial <- simGG(K = 400, pairs = 200, noff = 4, g = 15,
    nimm = 40, nimmG = seq(2, g-1, 1),		    # nimmG default value
    VAf = 1, VAi = 1, VRf = 1, VRi = 1,		    # all default values
    mup = 20, muf = 0, mui = 3, murf = 0, muri = 0, # mup and mui non-default values
    d_bvf = 0, d_bvi = 0, d_rf = 0, d_ri = 0)	    # all default values
 }
}

\keyword{datasets}

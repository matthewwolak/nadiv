\name{Mrode3}
\alias{Mrode3}
\docType{data}
\title{
Pedigree, from chapter 3 of Mrode (2005)}
\description{
An example pedigree with genetic groups and a data column
}
\usage{data(Mrode3)}
\format{
  A data frame with 10 observations on the following 8 variables.
  \describe{
    \item{\code{calf} }{a factor with levels indicating the unique genetic groups and individuals}
    \item{\code{dam} }{a numeric vector of maternal identities}
    \item{\code{sire} }{a numeric vector of paternal identities}
    \item{\code{damGG} }{a factor of maternal identities with genetic groups inserted instead of \code{NA}}
    \item{\code{sireGG} }{a factor of paternal identities with genetic groups inserted instead of \code{NA}}
    \item{\code{sex} }{a factor with levels \code{female} \code{male}}
    \item{\code{WWG} }{a numeric vector of pre-weaning weight gain (kg) for five beef calves}
  }
}
\source{Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding Values, 2nd ed.  Cambridge, MA: CABI Publishing.
}
\keyword{datasets}

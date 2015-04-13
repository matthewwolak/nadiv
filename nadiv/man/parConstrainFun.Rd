\name{parConstrainFun}
\alias{parConstrainFun}
\title{
Function used in the \code{proLik} function to produce a profile likelihood for a variance component}

\description{
Given a model object from \code{asreml} and a range of estimates of the parameter, the function will supply the likelihood ratio test statistic for the comparison of the full model to one where the parameter of interest is constrained.}

\usage{
parConstrainFun(x, parameters, full, fm2, comp, G)
}

\arguments{
  \item{x }{section of all parameter values to analyze}
  \item{parameters }{a value for which the log Likelihood of a model is to be calculated}
  \item{full }{the full model \code{asreml} object}
  \item{fm2 }{starting values for the full model}
  \item{comp }{which variance component to constrain}
  \item{G }{logical indicating if the component is part of the G structure}
}

\details{
Used internally in the \code{proLik} function to call \code{constrainFun}}

\author{\email{matthewwolak@gmail.com}}

\seealso{See Also \code{\link{proLik}}, \code{\link{constrainFun}}}

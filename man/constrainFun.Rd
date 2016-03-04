\name{constrainFun}
\alias{constrainFun}
\title{
Function used in conjunction with others to produce a profile likelihood for a variance component}

\description{
Given a model object from \code{asreml} and a range of estimates of the parameter, the function will supply the likelihood ratio test statistic for the comparison of the full model to one where the parameter of interest is constrained.}

\usage{
constrainFun(parameter.val, full, fm2, comp, G, mit = 600)
}

\arguments{
  \item{parameter.val }{a value for which the log-Likelihood of a model is to be calculated}
  \item{full }{the full model \code{asreml} object}
  \item{fm2 }{starting values for the full model}
  \item{comp }{which variance component to constrain}
  \item{G }{logical, indicating if the component is part of the G structure}
  \item{mit }{numeric, indicating maximum number of iterations for the constrained asreml model}
}

\details{
Used internally in the \code{proLik} function}

\author{\email{matthewwolak@gmail.com}}

\seealso{See also \code{\link{proLik}}}

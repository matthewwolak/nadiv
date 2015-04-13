\name{plot.proLik}
\alias{plot.proLik}
\alias{plot}

\title{Plot profile Likelihood}

\description{This function graphically depicts the output from a \code{proLik} object.}

\usage{
  \method{plot}{proLik}(x, CL = TRUE, alpha = 0.05, type = "l", ...)
}

\arguments{
  \item{x }{the output from \code{proLik}}
  \item{CL }{a logical indicating whether a line representing the Confidence Limit is to be drawn}
  \item{alpha }{if \code{CL} is TRUE, the confidence level at which \code{CL} is to be drawn}
  \item{type }{the type of plot to be generated, see arguments to \code{plot}}
  \item{\dots }{other arguments to \code{plot}}
}

\author{\email{matthewwolak@gmail.com}}


\seealso{\code{proLik}}

\examples{##See examples in the help file for proLik()}

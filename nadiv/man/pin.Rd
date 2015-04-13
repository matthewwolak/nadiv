\name{pin}
\alias{pin}
\title{
Approximate standard errors for linear functions of variance components
}
\description{
This function is similar to the pin calculations performed by the standalone ASReml.  This function, written by Ian White, applies the delta method for the estimation of approximate standard errors on linear functions of variance components from a REML mixed model
}
\usage{
pin(object, transform)
}
\arguments{
  \item{object }{
A list with at least the following elements: \code{gammas}, \code{gammas.type}, and \code{ai} from a REML mixed model
}
  \item{transform }{
A formula specifying the linear transformation of variance components to conduct
}
}
\details{
Object is intended to be an asreml-R model output.

The formula can use \code{V1,..., Vn} to specify any one of the \code{n} variance components.  These should be in the same order as they are in the object (e.g., see the row order of \code{summary(object)$varcomp} for asreml-R models.
}
\value{
A dataframe with row names corresponding to the operator on the left hand side of the \code{transform} formula and the entries corresponding to the \code{Estimate} and approximate \code{SE} of the linear transformation.
}

\author{
Ian White
}


\seealso{
See Also as \code{\link{aiCI}}, \code{\link{aiFun}}
}

\examples{
# Below is the heritability calculation for tait1 of the warcolak dataset
# Re-create the output from a basic, univariate animal model in asreml-R
   asrMod <- list(gammas = c(0.6383285, 1.00),
	gammas.type = c(2, 1),
	ai = c(0.0044461106, -0.0011754947, 0.0004424668))
   namevec <- c("ped(ID)!ped", "R!variance")
   names(asrMod[[1]]) <- names(asrMod[[2]]) <- namevec

   nadiv:::pin(asrMod, h2 ~ V1 / (V1 + V2))
}


\name{simGG}
\alias{simGG}
\title{
Genetic group pedigree and data simulation
}
\description{
Simulates a pedigree and phenotype for a focal population receiving immigrants. Genetic and environmental differences can be specified between the focal and immigrant populations. Further, these differences can have temporal trends. 
}

\usage{
simGG(K, pairs, noff, g,
	nimm = 2, nimmG = seq(2, g - 1, 1),
	VAf = 1, VAi = 1, VRf = 1, VRi = 1,
	mup = 20, muf = 0, mui = 0, murf = 0, muri = 0,
	d_bvf = 0, d_bvi = 0, d_rf = 0, d_ri = 0)
}

\arguments{
  \item{K }{
Integer number of individuals per generation, or the focal population carrying capacity
}
  \item{pairs }{
Integer number of mating pairs created by sampling with replacement from adults of a given generation
}
  \item{noff }{
Integer number of offspring each pair contributes to the next generation
}
  \item{g }{
Integer number of (non-overlapping) generations to simulate
}
  \item{nimm }{
Integer number of immigrants added to the population each generation of migration
}
  \item{nimmG }{
Sequence of integers for the generations in which immigrants arrive in the focal population
}
  \item{VAf }{
Numeric value for the expected additive genetic variance in the first generation of the focal population
}
  \item{VAi }{
Numeric value for the expected additive genetic variance in each generation of immigrants
}
  \item{VRf }{
Numeric value for the expected residual variance in the focal population
}
  \item{VRi }{
Numeric value for the expected residual variance in each generation of the immigrants
}
  \item{mup }{
Numeric value for the expected mean phenotypic value in the first generation of the founder population
}
  \item{muf }{
Numeric value for the expected mean breeding value in the first generation of the focal population
}
  \item{mui }{
Numeric value for the expected mean breeding value for the immigrants
}
  \item{murf}{
Numeric value for the expected mean residual (environmental) deviation in the first generation of the focal population
}
  \item{muri}{
Numeric value for the expected mean residual (environmental) deviation for the immigrants
}
  \item{d_bvf }{
Numeric value for the expected change between generations in the mean breeding value of the focal population. Sets the rate of genetic selection occurring across generations
}
  \item{d_bvi }{
Numeric value for the expected change between generations in the mean breeding value of the immigrant population each generation
}
  \item{d_rf }{
Numeric value for the expected change between generations in the mean residual (environmental) deviation of the focal population each generation
}
  \item{d_ri }{
Numeric value for the expected change between generations in the mean residual (environmental) deviation of the immigrant population each generation
}
}


\details{
Offspring total additive genetic values \code{u} are the average of their parents \code{u} plus a Mendelian sampling deviation drawn from a normal distribution with mean of 0 and variance equal to \eqn{0.5V_{A} (1 - f_{sd})} where \eqn{V_A} is \code{VAf} and \eqn{f_{sd}} is the average of the parents' coefficient of inbreeding \emph{f} (p. 447 Verrier et al. 1993). Each \sQuote{immigrant} (individual with unknown parents in generations >1) is given a total additive genetic effect that is drawn from a normal distribution with mean of \code{mui} and variance equal to \code{VAi}. Residual deviations are sampled for \sQuote{resident} and \sQuote{immigrant} populations separately, using normal distributions with means of \code{murf} and \code{muri}, respectively, and variances of \code{VRf} and \code{VRi}, respectively. Phenotypes are the sum of total additive genetic effects and residual deviations plus an overall mean \code{mup}.

Trends in total additive genetic effects and/or residual deviations can be specified for both the resident and immigrant populations. Trends in total additive genetic effects occurring in the immigrants, in the residual deviations occurring in the residents, and in the residual deviations occurring in the immigrants are produced by altering the mean each generation for the separate distribution from which these effects are each drawn. The change in mean over a generation is specified in units of standard deviations of the respective distributions (e.g., square roots of \code{VAi}, \code{VRf}, and \code{VRi}) and is set with \code{d_bvi}, \code{d_rf}, or \code{d_ri}, respectively.

Trends in total additive genetic effects for the resident population are produced by selecting individuals to be parents of the next generation according to their total additive genetic effects. Individuals are assigned probabilities of being selected as a parent of the next generation depending on how closely their total additive genetic effect matches an optimum value. The parameter \code{d_bvf} specifies how much the optimal total additive genetic effect changes per generation. \code{d_bvf} is multipled by the square root of \code{VAf} (i.e., the additive genetic standard deviation of the resident population) and the number of generations since selection began, and it is this quantity that is added to the base generation's mean total additive genetic effect (\code{muf}) to determine the optimal total additive genetic effect in a given generation. Individuals with total additive genetic effects closest to this optimum have a higher probability of being randomly sampled to be parents of the next generation. This represents selection directly on total additive genetic effects.
}


\value{
A \code{data.frame} with columns corresponding to:
  \item{id }{Integer for each individual's unique identifying code}
  \item{dam }{Integer indicating each individual's dam}
  \item{sire }{Integer indicating each individual's sire}
  \item{parAvgU }{Numeric value for the average of each individual's dam and sire additive genetic effects}
  \item{mendel }{Numeric value for each individual's Mendelian sampling deviate from the mid-parental total additive genetic value}
  \item{u }{Numeric value of each individual's total additive genetic effect}
  \item{r }{Numeric value of each individual's residual (environmental) deviation}
  \item{p }{Numeric value of each individual's phenotypic value}
  \item{is }{Integer of either \code{0} if an individual was born in the focal population or \code{1} if they were born in an immigrant population}
  \item{gen }{Integer value of the generation in which each individual was born}
}

\references{
Verrier, V., J.J. Colleau, and J.L. Foulley. 1993. Long-term effects of selection based on the animal model BLUP in a finite population. Theoretical and Applied Genetics. 87:446-454.
}

\author{\email{matthewwolak@gmail.com}}



\seealso{
  \code{\link{ggTutorial}}
}

\examples{
  \dontrun{
  # The dataset 'ggTutorial' was simulated as:
  set.seed(102)      		# seed used to simulate ggTutorial
  ggTutorial <- simGG(K = 400, pairs = 200, noff = 4, g = 15,
	nimm = 40,
	muf = 0, mui = 3)
  }
}



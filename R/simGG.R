#' Genetic group pedigree and data simulation
#' 
#' Simulates a pedigree and phenotype for a focal population receiving
#' immigrants. Genetic and environmental differences can be specified between
#' the focal and immigrant populations. Further, these differences can have
#' temporal trends.
#' 
#' Offspring total additive genetic values \code{u} are the average of their
#' parents \code{u} plus a Mendelian sampling deviation drawn from a normal
#' distribution with mean of 0 and variance equal to \eqn{0.5V_{A} (1 -
#' f_{sd})}{0.5V[A] (1 - f[sd])} where \eqn{V_A}{V[A]} is \code{VAf} and
#' \eqn{f_{sd}}{f[sd]} is the average of the parents' coefficient of inbreeding
#' \emph{f} (p. 447 Verrier et al. 1993). Each \sQuote{immigrant} (individual
#' with unknown parents in generations >1) is given a total additive genetic
#' effect that is drawn from a normal distribution with mean of \code{mui} and
#' variance equal to \code{VAi}. Residual deviations are sampled for
#' \sQuote{focal} and \sQuote{immigrant} populations separately, using normal
#' distributions with means of \code{murf} and \code{muri}, respectively, and
#' variances of \code{VRf} and \code{VRi}, respectively. Phenotypes are the sum
#' of total additive genetic effects and residual deviations plus an overall
#' mean \code{mup}.
#' 
#' Trends in total additive genetic effects and/or residual deviations can be
#' specified for both the focal and immigrant populations. Trends in total
#' additive genetic effects occurring in the immigrants, in the residual
#' deviations occurring in the focal population, and in the residual deviations
#' occurring in the immigrants are produced by altering the mean each
#' generation for the separate distribution from which these effects are each
#' drawn. The change in mean over a generation is specified in units of
#' standard deviations of the respective distributions (e.g., square roots of
#' \code{VAi}, \code{VRf}, and \code{VRi}) and is set with \code{d_bvi},
#' \code{d_rf}, or \code{d_ri}, respectively.
#' 
#' Trends in total additive genetic effects for the focal population are
#' produced by selecting individuals to be parents of the next generation
#' according to their \emph{predicted} total additive genetic effects.
#' Individuals are assigned probabilities of being selected as a parent of the
#' next generation depending on how closely their predicted total additive
#' genetic effect matches an optimum value. Probabilities are assigned:
#' \deqn{exp((\frac{-1}{2\sigma_{x}}) (x - \theta)^{2})}{exp((-1/(2 * sd(x))) *
#' (x - theta)^2)} where \code{x} is the vector of predicted total additive
#' genetic effects (\code{u}), \eqn{\sigma_{x}}{sd(x)} is the standard
#' deviation of \code{x}, and \eqn{\theta}{theta} is the optimum value.
#' Sampling is conducted with replacement based on these probabilities.
#' 
#' The parameter \code{d_bvf} specifies how much the optimal total additive
#' genetic effect changes per generation. The optimal total additive genetic
#' effect in a given generation is calculated as: \code{muf + d_bvf
#' *}\code{sqrt(VAf) * (i-2)}. Individuals with predicted total additive
#' genetic effects closest to this optimum have a higher probability of being
#' randomly sampled to be parents of the next generation. This represents
#' selection directly on predicted total additive genetic effects.
#' 
#' Total additive genetic effects are predicted for the first generation of
#' focal individuals and all immigrants using equation 1.3 in Mrode (2005,
#' p.3): \eqn{h^{2} * (phenotype_{i} - mean population phenotype)}. The
#' heritability is either \code{VAf} / (\code{VAf + VRf}) or \code{VAi} /
#' (\code{VAi + VRi}). Total additive genetic effects are predicted for all
#' other individuals using equation 1.9 in Mrode (2005, p. 10) - or as the
#' average of each individual's parents' predicted total additive genetic
#' effects.
#' 
#' @param K Integer number of individuals per generation, or the focal
#'   population carrying capacity
#' @param pairs Integer number of mating pairs created by sampling with
#'   replacement from adults of a given generation
#' @param noff Integer number of offspring each pair contributes to the next
#'   generation
#' @param g Integer number of (non-overlapping) generations to simulate
#' @param nimm Integer number of immigrants added to the population each
#'   generation of migration
#' @param nimmG Sequence of integers for the generations in which immigrants
#'   arrive in the focal population
#' @param VAf Numeric value for the expected additive genetic variance in the
#'   first generation of the focal population - the founders
#' @param VAi Numeric value for the expected additive genetic variance in each
#'   generation of immigrants
#' @param VRf Numeric value for the expected residual variance in the focal
#'   population
#' @param VRi Numeric value for the expected residual variance in each
#'   generation of the immigrants
#' @param mup Numeric value for the expected mean phenotypic value in the first
#'   generation of the focal population - the founders
#' @param muf Numeric value for the expected mean breeding value in the first
#'   generation of the focal population - the founders
#' @param mui Numeric value for the expected mean breeding value for the
#'   immigrants
#' @param murf Numeric value for the expected mean residual (environmental)
#'   deviation in the first generation of the focal population - the founders
#' @param muri Numeric value for the expected mean residual (environmental)
#'   deviation for the immigrants
#' @param d_bvf Numeric value for the expected change between generations in
#'   the mean breeding value of the focal population. Sets the rate of genetic
#'   selection occurring across generations
#' @param d_bvi Numeric value for the expected change between generations in
#'   the mean breeding value of the immigrant population each generation
#' @param d_rf Numeric value for the expected change between generations in the
#'   mean residual (environmental) deviation of the focal population each
#'   generation
#' @param d_ri Numeric value for the expected change between generations in the
#'   mean residual (environmental) deviation of the immigrant population each
#'   generation
#'
#' @return A \code{data.frame} with columns corresponding to:
#'   \describe{
#'     \item{id }{Integer for each individual's unique identifying code}
#'     \item{dam }{Integer indicating each individual's dam}
#'     \item{sire }{Integer indicating each individual's sire}
#'     \item{parAvgU }{Numeric value for the average of each individual's dam 
#'       and sire additive genetic effects}
#'     \item{mendel }{Numeric value for each individual's Mendelian sampling 
#'       deviate from the mid-parental total additive genetic value}
#'     \item{u }{Numeric value of each individual's total additive genetic 
#'       effect}
#'     \item{r }{Numeric value of each individual's residual (environmental) 
#'       deviation}
#'     \item{p }{Numeric value of each individual's phenotypic value}
#'     \item{pred.u }{Numeric value of each individual's predicted total 
#'       additive genetic effect}
#'     \item{is }{Integer of either \code{0} if an individual was born in the 
#'       focal population or \code{1} if they were born in an immigrant 
#'       population}
#'     \item{gen }{Integer value of the generation in which each individual was 
#'       born}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{ggTutorial}}
#' @references Verrier, V., J.J. Colleau, and J.L. Foulley. 1993. Long-term
#' effects of selection based on the animal model BLUP in a finite population.
#' Theoretical and Applied Genetics. 87:446-454.
#' 
#' Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding
#' Values, 2nd ed.  Cambridge, MA: CABI Publishing.
#' @examples
#' 
#'   \dontrun{
#'   # The dataset 'ggTutorial' was simulated as:
#'   set.seed(102)      		# seed used to simulate ggTutorial
#'   ggTutorial <- simGG(K = 400, pairs = 200, noff = 4, g = 15,
#' 	nimm = 40,
#' 	muf = 0, mui = 3)
#'   }
#' 
#' # Use genetic group methods to approximate the breeding values for ggTutorial
#'   ## First, construct a pedigree with genetic groups
#'   ggPed <- ggTutorial[, c("id", "dam", "sire", "is", "gen")]
#'   naPar <- which(is.na(ggPed[, 2]))
#'   ggPed$GG <- rep("NA", nrow(ggPed))
#'     # 'focal' population genetic group = "foc0" and 'immigrant' = "g1"
#'     # obtained by pasting "foc" & "g" with immigrant status "0" or "1", respectively
#'     ggPed$GG[naPar] <- as.character(ggPed$is[naPar])
#'     ggPed$GG[ggPed$GG == "0"] <- paste0("foc", ggPed$GG[ggPed$GG == "0"])
#'     ggPed$GG[ggPed$GG == "1"] <- paste0("g", ggPed$GG[ggPed$GG == "1"])
#'   ggPed[naPar, 2:3] <- ggPed[naPar, "GG"]
#' 
#'   ## Now create the Q matrix
#'   Q <- ggcontrib(ggPed[, 1:3], ggroups = c("foc0", "g1"))
#' 
#'   ## obtain the true values of the genetic group means
#'   foc0_mean <- mean(ggTutorial$u[which(ggTutorial$gen == 1 & ggTutorial$is == 0)])
#'   g1_mean <- mean(ggTutorial$u[which(ggTutorial$is == 1)])
#'   g_exp <- matrix(c(foc0_mean, g1_mean), ncol = 1)
#' 
#'   ## breeding values (a) are:
#'   ### tot. add. gen. effects (u) minus genetic group effects for each individual (Qg):
#'   a <- ggTutorial$u - Q %*% g_exp
#' 
#' 
#' @export
simGG2 <- function(K, pairs, noff, g,
	nimm = 2, nimmG = seq(2, g-1, 1),
	VAf = 1, VAi = 1, VRf = 1, VRi = 1,
	mup = 20, muf = 0, mui = 0, murf = 0, muri = 0,
	d_bvf = 0, d_bvi = 0, d_rf = 0, d_ri = 0){

  if(pairs*2 > K) stop("pairs must be less than half of K")
  if(nimmG[1] == 1) stop("immigrants cannot arrive in the first generation")
  N <- pairs*noff
  da <- array(NA, dim = c(K, 12, g))
  dimnames(da) <- list(NULL,
     c("id", "dam", "sire", "parAvgU", "mendel", "u", "r", "p", "pred.u", "is", "gen", "q_imm"), seq(g))
  da[, "id", ] <- seq(K*g)
  # Assume last nimm rows in each generation are the immigrants: 
  ## 1=immigrant & 0=NOT immigrant
  ## initially only migrants have non-zero immigrant genetic group contribution (`q_imm`)
  da[, c("is", "q_imm"), ] <- 0
  da[(K-nimm+1):K, c("is", "q_imm"), nimmG] <- 1
  da[, "gen", ] <- rep(seq(g), each = K)
  # Create standard normals for generation 1
  da[, "u", 1] <- rnorm(K, muf, sqrt(VAf))
  da[, "r", 1] <- rnorm(K, murf, sqrt(VRf))
  # Calculate phenotypes for generation 1
  da[, "p", 1] <- mup + rowSums(da[, c("u", "r"), 1])
  # Calculate predicted total additive genetic effects for generation 1
  ## equation 1.3 in Mrode (2005, p. 3)
  da[, "pred.u", 1] <- (VAf / (VAf + VRf)) * (da[, "p", 1] - mean(da[, "p", 1]))

  ########
  # Mating
  for(i in 2:g){
    if(i %in% nimmG) Knimm <- K-nimm else Knimm <- K
    # Parents ranked: how close pred.u is to the next generation's mean u
    ## Define function to do this
    prFun <- function(x, theta = 0){
      exp((-1 / (2*sd(x))) * (x - theta)^2)
    }
    ### End function definition
    # Arbitrarily choose sexes (for assignment as female=0 or male=1 parent)
    sexvec <- vector("integer", length = K)
    sexvec[sample.int(K, size = K*0.5, replace = FALSE)] <- 1
    if(d_bvf == 0){
      prin0 <- prin1 <- NULL
    } else{
        prin0 <- prFun(da[sexvec == 0, "pred.u", i-1],
	  muf + d_bvf*sqrt(VAf)*(i-2)) # no selection in first generation hence i-2
        prin1 <- prFun(da[sexvec == 1, "pred.u", i-1],
	  muf + d_bvf*sqrt(VAf)*(i-2)) # no selection in first generation hence i-2
      }
    # Sampling WITH replacement to assign parents
    pool0 <- sample(x = da[sexvec == 0, "id", i-1], size = pairs, replace = TRUE,
	prob = prin0)
    pool1 <- sample(x = da[sexvec == 1, "id", i-1], size = pairs, replace = TRUE,
	prob = prin1)
    iOff <- matrix(rep(c(pool0, pool1), each = noff), ncol = 2)
    da[1:Knimm, c("dam", "sire"), i] <- iOff[sort(sample(seq(nrow(iOff)),
	size = Knimm, replace = FALSE)), ]
    # Average of parent total additive genetic effects
    da[1:Knimm, "parAvgU", i] <- rowMeans(matrix(da[match(da[1:Knimm, c("dam", "sire"), i],
	da[, "id", i-1]), "u", i-1], ncol = 2, byrow = FALSE))
    # Average of parent immigrant genetic group contributions (`q_imm`)
    da[1:Knimm, "q_imm", i] <- rowMeans(matrix(da[match(da[1:Knimm, c("dam", "sire"), i],
	da[, "id", i-1]), "q_imm", i-1], ncol = 2, byrow = FALSE))

    # Assign Mendelian sampling variation
    # Within-family additive genetic variance
    ## p. 447, second eqn. in Verrier, Colleau, & Foulley. 1993. Theor. Appl. Genetics
    # Genetic group specific VAs, so group-specific sampling variances
    ## Splitting breeding values, but don't have to do that 
    ### Just split Mendelian sampling deviations by Muff et al. eqn. (10)
    #### Muff et al. eqn. (10) gives proportion of group VA for distribution of
    ##### group-specific Mendelian sampling deviations
    if(i == 2){  #<-- all ancestors from same genetic group: no need to split BVs
      da[1:Knimm, "mendel", i] <- rnorm(Knimm, 0, sqrt(0.5 * VAf))
    } else{
      tmpPed <- prunePed(data.frame(apply(da[, 1:3, 1:i], MARGIN = 2,
		FUN = function(x){x})), as.character(unique(c(da[1:Knimm, "id", i]))))
      igen_dii <- makeAinv(tmpPed)$dii[match(as.character(c(da[1:Knimm, "id", i])),
				tmpPed[, 1])]
#TODO CHECK this      # If homogeneous genetic group VAs, then don't need to split
      if(VAf == VAi){
        da[1:Knimm, "mendel", i] <- rnorm(Knimm, 0, sqrt(igen_dii * VAf))
      } else{
        ## Split Mendelian sampling deviations (Muff et al. eqn. 10)
        ### Total Mendelian sampling deviation is:
        #### rnorm(1, 0, sqrt(d_ii^(f) * VAf)) + rnorm(1, 0, sqrt(d_ii^(i) * VAi))
          da[1:Knimm, "mendel", i] <- rnorm(Knimm, 0,
		sqrt((1 - (1-da[1:Knimm, "q_imm", i]) * (1 - igen_dii)) * VAf)) +
	    rnorm(Knimm, 0, sqrt((1 - da[1:Knimm, "q_imm", i] * (1 - igen_dii)) * VAi))
        }
      }
    # Total additive genetic effects
    da[1:Knimm, "u", i] <- rowSums(da[1:Knimm, c("parAvgU", "mendel"), i])
    # Calculate predicted total additive genetic effects
    ## average of parents' predicted total additive genetic effects
    ## equation 1.9 in Mrode (2005, p. 10)
    da[1:Knimm, "pred.u", i] <- rowMeans(matrix(da[match(da[1:Knimm, c("dam", "sire"), i],
	da[, "id", i-1]), "pred.u", i-1], ncol = 2, byrow = FALSE))
    # Residual deviations
    da[1:Knimm, "r", i] <- rnorm(Knimm, murf + d_rf*sqrt(VRf)*(i-1), sqrt(VRf))
    #############
    # Immigrants: assumed outbred
    if(i %in% nimmG){
      # no trend in first generation (hence i-2)
      da[(Knimm+1):K, "u", i] <- rnorm(nimm, mui + d_bvi*sqrt(VAi)*(i-2), sqrt(VAi))
      # no trend in first generation (hence i-2)
      da[(Knimm+1):K, "r", i] <- rnorm(nimm, muri + d_ri*sqrt(VRi)*(i-2), sqrt(VRi))
    }
    #############
    # All individuals in generation 'i'
    # Calculate phenotypes
    da[, "p", i] <- mup + rowSums(da[, c("u", "r"), i])
    # Calculate predicted total additive genetic effects for immigrants
    ## equation 1.3 in Mrode (2005, p. 3)
    if(i %in% nimmG){
      da[(Knimm+1):K, "pred.u", i] <- (VAi / (VAi + VRi)) * (da[(Knimm+1):K, "p", i] - mean(da[(Knimm+1):K, "p", i]))
    }    
  } 
  # create a data.frame out of the array
  df <- data.frame(apply(da, MARGIN = 2, FUN = function(x){x}))
 df
} 

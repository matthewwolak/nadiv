#' Pedigree adapted from Fikse 2009 with genetic groups and fuzzy classification
#' 
#' @format A \code{data.frame} with 16 observations on the following 11 variables:
#'   \describe{
#'     \item{id }{a factor with levels indicating the unique individuals 
#'       (including phantom parents) and genetic groups}
#'     \item{dam }{a factor of observed maternal identities}
#'     \item{sire }{a factor vector of observed paternal identities}
#'     \item{damGG }{a factor of maternal identities with genetic groups
#'       inserted instead of \code{NA}}
#'     \item{sireGG }{a factor of paternal identities with genetic groups
#'       inserted instead of \code{NA}}
#'     \item{phantomDam }{a factor of maternal identities with phantom parents
#'       inserted instead of \code{NA}}
#'     \item{phantomSire }{a factor of paternal identities with phantom parents
#'       inserted instead of \code{NA}}
#'     \item{group }{a factor of genetic groups to which each phantom parent
#'       belongs}
#'     \item{g1 }{a numeric vector with probabilities of group \code{g1}
#'       membership for each phantom parent}
#'     \item{g2 }{a numeric vector with probabilities of group \code{g2}
#'       membership for each phantom parent}
#'     \item{g3 }{a numeric vector with probabilities of group \code{g3}
#'       membership for each phantom parent}
#'   }
#'
#' @docType data
#' @source Fikse, F. 2009. Fuzzy classification of phantom parent groups in an
#'   animal model. Genetics Selection Evolution 41:42.
#' @keywords datasets
#' @examples
#'   data(F2009)
#'   str(F2009)
"F2009"




#' Pedigree, adapted from Table 1 in Fernando & Grossman (1990)
#' 
#' @format A \code{data.frame} with 8 observations on the following 4 variables:
#'   \describe{
#'     \item{id }{a factor with levels \code{1} \code{2} \code{3} \code{4} 
#'       \code{5} \code{6} \code{7} \code{8}}
#'     \item{dam }{a factor with levels \code{2} \code{4} \code{6}}
#'     \item{sire }{a factor with levels \code{1} \code{3} \code{5}}
#'     \item{sex }{a factor with levels \code{0} \code{1}}
#'   }
#'
#' @docType data
#' @source Fernando, R.L. & M. Grossman. 1990. Genetic evaluation with
#' autosomal and X-chromosomal inheritance. Theoretical and Applied Genetics
#' 80:75-80.
#' @keywords datasets
#' @examples
#'   data(FG90)
#'   str(FG90)
"FG90"





#' Simulated dataset used to analyze data with genetic group animal models
#' 
#' The dataset was simulated using the \code{\link{simGG}} function so that the
#' pedigree contains a base population comprised of founders and non-founder
#' immigrants. These data are then used in the main manuscript and tutorials
#' accompanying Wolak & Reid (2017).
#' 
#' The dataset was simulated as described in the \sQuote{examples} section
#' using the \code{\link{simGG}} function. Full details of the function and
#' dataset can be found in Wolak & Reid (2017).
#' 
#' The \code{data.frame} contains 6000 individuals across 15 generations. In
#' each generation, the carrying capacity is limited to 400 individuals, the
#' number of mating pairs limited to 200 pairs, and 40 immigrants per
#' generation arrive starting in the second generation.
#' 
#' The breeding values of the founders are drawn from a normal distribution
#' with an expected mean of 0 and a variance of 1. The breeding values of all
#' immigrants are drawn from a normal distribution with an expected mean of 3
#' and variance of 1. Consequently, the expected difference between mean
#' breeding values in the founders and immigrants is 3. All individuals are
#' assigned a residual (environmental) deviation that is drawn from a normal
#' distribution with an expected mean of 0 and variance of 1.
#' 
#' @format A \code{data.frame} with 6000 observations on the following 10 
#'   variables:
#'   \describe{
#'     \item{id }{an integer vector specifying the 6000 unique individual
#'       identities}
#'     \item{dam }{an integer vector specifying the unique dam for each 
#'       individual}
#'     \item{sire }{an integer vector specifying the unique sire for each 
#'       individual}
#'     \item{parAvgU }{a numeric vector of the average autosomal total additive 
#'       genetic effects (\code{u}) of each individual's parents}
#'     \item{mendel }{a numeric vector of the Mendelian sampling deviations 
#'       from \code{parAvgU} autosomal total additive genetic effects that is 
#'       unique to each individual}
#'     \item{u }{a numeric vector of the total autosomal additive genetic 
#'       effects underlying \code{p}}
#'     \item{r }{a numeric vector of the residual (environmental) effects 
#'       underlying \code{p}}
#'     \item{p }{a numeric vector of phenotypic values}
#'     \item{is}{an integer vector with \code{0} for individuals born in the 
#'       focal population and \code{1} for individuals born outside of the 
#'       focal population, but immigrated}
#'     \item{gen }{an integer vector specifying the generation in which each 
#'       individual was born}
#'   }
#'
#' @docType data
#' @source Wolak, M.E. & J.M. 2017. Accounting for genetic differences among
#'   unknown parents in microevolutionary studies: how to include genetic
#'   groups in quantitative genetic animal models. Journal of Animal Ecology
#'   86:7-20. doi:10.1111/1365-2656.12597
#' @keywords datasets
#' @examples
#' 
#'  \dontrun{
#'   rm(list = ls())
#'   set.seed(102)     #<-- seed value used originally
#'   library(nadiv)
#'   # create data using `simGG()`
#'   ggTutorial <- simGG(K = 400, pairs = 200, noff = 4, g = 15,
#'     nimm = 40, nimmG = seq(2, g-1, 1),		    # nimmG default value
#'     VAf = 1, VAi = 1, VRf = 1, VRi = 1,		    # all default values
#'     mup = 20, muf = 0, mui = 3, murf = 0, muri = 0, # mup and mui non-default values
#'     d_bvf = 0, d_bvi = 0, d_rf = 0, d_ri = 0)	    # all default values
#'  }
#' 
"ggTutorial"





#' Pedigree from Table 2.1 of Mrode (2005)
#' 
#' @format A \code{data.frame} with 6 observations on the following 3 variables:
#'   \describe{
#'     \item{id }{a numeric vector}
#'     \item{dam }{a numeric vector}
#'     \item{sire }{a numeric vector}
#'   }
#'
#' @docType data
#' @source Mrode, R.A. 2005. Linear Models for the Prediction of Animal
#' Breeding Values, 2nd ed.  Cambridge, MA: CABI Publishing.
#' @keywords datasets
#' @examples
#    data(Mrode2)
#'   str(Mrode2)
"Mrode2"





#' Pedigree, from chapter 3 of Mrode (2005) with genetic groups and a trait column
#' 
#' @format A \code{data.frame} with 10 observations on the following 8 variables:
#'   \describe{
#'     \item{calf }{a factor with levels indicating the unique genetic groups 
#'       and individuals}
#'     \item{dam }{a numeric vector of maternal identities}
#'     \item{sire }{a numeric vector of paternal identities}
#'     \item{damGG }{a factor of maternal identities with genetic groups 
#'       inserted instead of \code{NA}}
#'     \item{sireGG }{a factor of paternal identities with genetic groups 
#'       inserted instead of \code{NA}}
#'     \item{sex }{a factor with levels \code{female} \code{male}}
#'     \item{WWG }{a numeric vector of pre-weaning weight gain (kg) for five 
#'       beef calves}
#'   }
#'
#' @docType data
#' @source Mrode, R.A. 2005. Linear Models for the Prediction of Animal
#' Breeding Values, 2nd ed.  Cambridge, MA: CABI Publishing.
#' @keywords datasets
#' @examples
#'   data(Mrode3)
#'   str(Mrode3)
"Mrode3"





#' Pedigree, adapted from example 9.1 of Mrode (2005)
#' 
#' @format A \code{data.frame} with 12 observations on the following 3 variables:
#'   \describe{
#'     \item{pig }{a numeric vector}
#'     \item{dam }{a numeric vector}
#'     \item{sire }{a numeric vector}
#'   }
#'
#' @docType data
#' @source Mrode, R.A. 2005. Linear Models for the Prediction of Animal
#' Breeding Values, 2nd ed.  Cambridge, MA: CABI Publishing.
#' @keywords datasets
#' @examples
#'   data(Mrode9)
#'   str(Mrode9)
"Mrode9"





#' Pedigree with genetic groups adapted from Quaas (1988) equation [5]
#' 
#' @format A \code{data.frame} with 11 observations on the following 8 variables:
#'   \describe{
#'     \item{id }{a factor with levels indicating the unique individuals 
#'       (including phantom parents) and genetic groups}
#'     \item{dam }{a factor of observed maternal identities}
#'     \item{sire }{a factor vector of observed paternal identities}
#'     \item{damGG }{a factor of maternal identities with genetic groups 
#'       inserted instead of \code{NA}}
#'     \item{sireGG }{a factor of paternal identities with genetic groups 
#'       inserted instead of \code{NA}}
#'     \item{phantomDam }{a factor of maternal identities with phantom parents 
#'       inserted instead of \code{NA}}
#'     \item{phantomSire }{a factor of paternal identities with phantom parents 
#'       inserted instead of \code{NA}}
#'     \item{group }{a factor of genetic groups to which each phantom parent
#'       belongs}
#'   }
#' 
#' @docType data
#' @source Quaas, R.L. 1988. Additive genetic model with groups and
#' relationships. Journal of Dairy Science 71:1338-1345.
#' @keywords datasets
#' @examples
#'   data(Q1988)
#'   str(Q1988)
"Q1988"





#' Pedigree and phenotypic values for a mythical population of Warcolaks
#' 
#' A two trait example pedigree from the three generation breeding design of
#' Fairbairn & Roff (2006) with two un-correlated traits.
#' 
#' Unique sets of relatives are specified for a three generation breeding
#' design (Fairbairn & Roff, 2006).  Each set contains 72 individuals. This
#' pedigree reflects an experiment which produces 75 of these basic sets from
#' Fairbairn & Roff's design. The pedigree was created using
#' \code{simPedDFC()}.
#' 
#' The dataset was simulated to have two un-correlated traits with different
#' genetic architectures (see \code{examples} below). The trait means are both
#' equal to 1 for males and 2 for females. The additive genetic, dominance
#' genetic, and environmental (or residual) variances for both \code{trait1}
#' and \code{trait2} are 0.4, 0.3, & 0.3, respectively. However, the additive
#' genetic variance for \code{trait2} can be further decomposed to autosomal
#' additive genetic variance (0.3) and X-linked additive genetic variance (0.1;
#' assuming the \sQuote{no global dosage compensation} mechanism).
#' 
#' Females and males have equal variances (except for sex-chromosomal additive
#' genetic variance, where by definition, males have half of the additive
#' genetic variance as females; Wolak 2013) and a between-sex correlation of
#' one for all genetic and residual effects (except the cross-sex residual
#' covariance=0). All random effects were drawn from multivariate random normal
#' distributions [e.g., autosomal additive genetic effects: N ~ (0,
#' kronecker(A, G))] with means of zero and (co)variances equal to the product
#' of the expected sex-specific (co)variances (e.g., G) and the relatedness (or
#' incidence) matrix (e.g., A).
#' 
#' The actual variance in random effects will vary slightly from the amount
#' specified in the simulation, because of Monte Carlo error. Thus, the random
#' effects have been included as separate columns in the dataset. See
#' \code{examples} below for the code that generated the dataset.
#'
#' @format A \code{data.frame} with 5400 observations on the following 13 variables:
#'   \describe{
#'     \item{ID }{a factor specifying 5400 unique individual identities}
#'     \item{Dam }{a factor specifying the unique Dam for each individual}
#'     \item{Sire }{a factor specifying the unique Sire for each individual}
#'     \item{sex }{a factor specifying \dQuote{M} if the individual is a male 
#'       and \dQuote{F} if it is a female}
#'     \item{trait1 }{a numeric vector of phenotypic values: see 
#'       \sQuote{Details}}
#'     \item{trait2 }{a numeric vector of phenotypic values: see 
#'       \sQuote{Details}}
#'     \item{t1_a }{a numeric vector of the autosomal additive genetic effects 
#'       underlying \sQuote{trait1}}
#'     \item{t2_a }{a numeric vector of the autosomal additive genetic effects 
#'       underlying \sQuote{trait2}}
#'     \item{t2_s }{a numeric vector of the sex-chromosomal additive genetic 
#'       effects underlying \sQuote{trait2}}
#'     \item{t1_d }{a numeric vector of the autosomal dominance genetic effects 
#'       underlying \sQuote{trait1}}
#'     \item{t2_d }{a numeric vector of the autosomal dominance genetic effects 
#'       underlying \sQuote{trait2}}
#'     \item{t2_r }{a numeric vector of the residual (environmental) effects 
#'       underlying \sQuote{trait1}}
#'     \item{t2_r }{a numeric vector of the residual (environmental) effects 
#'       underlying \sQuote{trait2}}
#'   }
#' @note Before nadiv version 2.14.0, the \code{warcolak} dataset used a 0/1
#' coding for \sQuote{sex} and did not contain the random effects.
#' 
#' @docType data
#' @references Fairbairn, D.J. & Roff, D.A. 2006. The quantitative genetics of
#' sexual dimorphism: assessing the importance of sex-linkage. Heredity 97,
#' 319-328.
#' 
#' Wolak, M.E. 2013. The Quantitative Genetics of Sexual Differences: New
#' Methodologies and an Empirical Investigation of Sex-Linked, Sex-Specific,
#' Non-Additive, and Epigenetic Effects. Ph.D. Dissertation. University of
#' California Riverside.
#' @keywords datasets
#' @examples
#' 
#'  \dontrun{
#'   rm(list = ls())
#'   set.seed(101)
#'   library(nadiv)
#'   # create pedigree
#'   warcolak <- simPedDFC(F = 75, gpn = 4, fsn = 4, s = 2)
#'   names(warcolak)[1:3] <- c("ID", "Dam", "Sire")
#'   warcolak$trait2 <- warcolak$trait1 <- NA
#' 
#'   # Define covariance matrices for random effects:
#'   ## Autosomal additive genetic (trait1)
#'   Ga_t1 <- matrix(c(0.4, rep(0.399999, 2), 0.4), 2, 2)
#'   ## Autosomal additive genetic (trait2)
#'   Ga_t2 <- matrix(c(0.3, rep(0.299999, 2), 0.3), 2, 2)
#'   ## Sex-chromosomal additive genetic (trait2)
#'   Gs_t2 <- matrix(c(0.1, rep(0.099999, 2), 0.1), 2, 2)
#'   ## Autosomal dominance genetic
#'   Gd <- matrix(c(0.3, rep(0.299999, 2), 0.3), 2, 2)
#'   ## Environmental (or residual)
#'   ### Assumes no correlated environmental effects between sexes
#'   R <- diag(c(0.3, 0.3))
#' 
#'   ## define variables to be re-used
#'   pedn <- nrow(warcolak)
#'   # Female (homogametic sex chromosomes) in first column
#'   # Male (heterogametic sex chromosomes) in second column
#'   sexcol <- as.integer(warcolak$sex)
#' 
#'   # Create random effects
#'   ## Additive genetic
#'   ### trait1 autosomal
#'   tmp_a <- grfx(pedn, G = Ga_t1, incidence = makeA(warcolak[, 1:3]))
#'     var(tmp_a)
#'   warcolak$t1_a <- tmp_a[cbind(seq(pedn), sexcol)]
#'   ### trait2 autosomal
#'   tmp_a <- grfx(pedn, G = Ga_t2, incidence = makeA(warcolak[, 1:3]))
#'     var(tmp_a)
#'   warcolak$t2_a <- tmp_a[cbind(seq(pedn), sexcol)]
#'   ### trait2 sex-chromosomal
#'   tmp_s <- grfx(pedn, G = Gs_t2, incidence = makeS(warcolak[, 1:4],
#' 	heterogametic = "M", DosageComp = "ngdc", returnS = TRUE)$S)
#'     matrix(c(var(tmp_s[which(sexcol == 1), 1]),
#' 	rep(cov(tmp_s), 2), var(tmp_s[which(sexcol == 2), 2])), 2, 2)
#'     # NOTE above should be: covar = male var = 0.5* female var 
#'   warcolak$t2_s <- tmp_s[cbind(seq(pedn), sexcol)]
#' 
#'   ## Dominance genetic
#'   ### trait1 
#'   tmp_d <- grfx(pedn, G = Gd, incidence = makeD(warcolak[, 1:3], invertD = FALSE)$D)
#'     var(tmp_d)
#'   warcolak$t1_d <- tmp_d[cbind(seq(pedn), sexcol)]
#'   ### trait2
#'   tmp_d <- grfx(pedn, G = Gd, incidence = makeD(warcolak[, 1:3], invertD = FALSE)$D)
#'     var(tmp_d)
#'   warcolak$t2_d <- tmp_d[cbind(seq(pedn), sexcol)]
#' 
#'   ## Residual
#'   ### trait1
#'   tmp_r <- grfx(pedn, G = R, incidence = NULL) # warning of identity matrix
#'     var(tmp_r)
#'   warcolak$t1_r <- tmp_r[cbind(seq(pedn), sexcol)]
#'   ### trait2
#'   tmp_r <- grfx(pedn, G = R, incidence = NULL) # warning of identity matrix
#'     var(tmp_r)
#'   warcolak$t2_r <- tmp_r[cbind(seq(pedn), sexcol)]
#' 
#'   # Sum random effects and add sex-specific means to get phenotypes
#'   ## females have slightly greater mean trait values
#'   warcolak$trait1 <- 1 + (-1*sexcol + 2) + rowSums(warcolak[, c("t1_a", "t1_d", "t1_r")])
#'   warcolak$trait2 <- 1 + (-1*sexcol + 2) + rowSums(warcolak[, c("t2_a", "t2_s", "t2_d", "t2_r")])
#'  }
#' 
"warcolak"




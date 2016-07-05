simGG <- function(K, pairs, noff, g,
	nimm = 2, nimmG = seq(2, g-1, 1),
	VAf = 1, VAi = 1, VRf = 1, VRi = 1,
	mup = 20, muf = 0, mui = 0, murf = 0, muri = 0,
	d_bvf = 0, d_bvi = 0, d_rf = 0, d_ri = 0){

  if(pairs*2 > K) stop("pairs must be less than half of K")
  if(nimmG[1] == 1) stop("immigrants cannot arrive in the first generation")
  N <- pairs*noff
  da <- array(NA, dim = c(K, 11, g))
  dimnames(da) <- list(NULL,
     c("id", "dam", "sire", "parAvgU", "mendel", "u", "r", "p", "pred.u", "is", "gen"), seq(g))
  da[, "id", ] <- seq(K*g)
  # Assume last nimm rows in each generation are the immigrants: 
  ## 1=immigrant & 0=NOT immigrant
  da[, "is", ] <- 0
  da[(K-nimm+1):K, "is", nimmG] <- 1
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
    # Assign Mendelian sampling variation
    # Within-family additive genetic variance
    ## p. 447, second eqn. in Verrier, Colleau, & Foulley. 1993. Theor. Appl. Genetics
    if(i == 2){
      da[1:Knimm, "mendel", i] <- rnorm(Knimm, 0, sqrt(0.5 * VAf))
    } else{
      # Average of parent inbreeding coefficients
        tmpPed <- prunePed(data.frame(apply(da[, 1:3, 1:i-1], MARGIN = 2,
		FUN = function(x){x})), as.character(unique(c(da[1:Knimm, c("dam", "sire"), i]))))
        tmpParF <- makeAinv(tmpPed)$f
        da[1:Knimm, "mendel", i ] <- rnorm(Knimm, 0,
		sqrt(0.5 * VAf * (1 - rowMeans(matrix(tmpParF[match(as.character(c(da[1:Knimm,
		c("dam", "sire"), i])), tmpPed[, 1])], ncol = 2, byrow = FALSE)))))
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

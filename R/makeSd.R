#rm(list = ls()); library(nadiv)


makeSd <- function(pedigree, heterogametic,
	DosageComp = c(NULL, "ngdc", "hori", "hedo", "hoha", "hopi"),
	parallel = FALSE, ncores = getOption("mc.cores", 2L),
	invertSd = TRUE, returnS = FALSE, det = TRUE){


  if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")
  nPed <- numPed(pedigree[, 1:3])

  dc.model <- match.arg(DosageComp)
  if(is.null(dc.model)){
    warning("Assuming 'ngdc' dosage compensation model")
    dc.model <- "ngdc"
  }
  #TODO check if can have dominance under hori
  #TODO check if any inbreeding effects in hets for hedo (don't think so)
  #TODO check if any inbreeding effects in hoha
  if(dc.model == "hopi"){
    #TODO check this is right
    cat("Assume that dominance allelic interactions cannot occur under 'hopi'\n")
    return(NULL)
  }
if(dc.model != "ngdc"){ #FIXME temporarily only allow ngdc for now
  stop("Currently only supporting 'ngdc' model")
}

  Sout <- makeS(cbind(nPed, pedigree[, 4]), heterogametic = heterogametic,
	DosageComp = dc.model, returnS = returnS)
  # makeA() returns `dsCMatrix`, but S is `dgCMatrix`
  ## makeD()-like code expects symmetric matrix
  S <- forceSymmetric(Sout$S)

  damsex <- pedigree[unique(nPed[, 2])[-1], 4]
  if(any(damsex == heterogametic)){
    pedname <- names(pedigree)
    pedigree <- pedigree[, c(1,3,2,4)]
    names(pedigree) <- pedname
    nPed <- numPed(pedigree[, 1:3])
  }

  sex <- rep(-998, dim(pedigree)[1L])
  sex[homs <- which(pedigree[,4] != heterogametic)] <- 1
  sex[hets <- which(pedigree[,4] == heterogametic)] <- 0
  N <- dim(nPed)[1L]
#FIXME turned off next check so can test parallel=TRUE on small pedigrees
#  if(parallel){
#    if(length(S@x)/ncores < 10){
#      warning("pedigree too small - 'parallel' set to FALSE instead")
#      parallel <- FALSE
#    }
#  }

  if(!parallel){
     cat("starting to make Sd...")

     Cout <- .C("sdij",
                as.integer(nPed[, 2] - 1), 
		as.integer(nPed[, 3] - 1), 
		as.integer(S@i), 			
		as.integer(S@p),                        
		as.double(S@x/2),                       
		as.integer(N),                           
		as.double(rep(0, length(S@x))),         
		as.integer(rep(0, length(S@i))),        
                as.integer(rep(0, N)),                  
		as.integer(0),
		as.integer(sex))	                        

     Sd <- Matrix(0, N, N,
	sparse = TRUE, dimnames = list(as.character(pedigree[, 1]), NULL))
     Sd@uplo <- "U"
     Sd@i <- Cout[[8]][1:Cout[[10]]]
     Sd@p <- c(Cout[[9]], Cout[[10]])
     Sd@x <- Cout[[7]][1:Cout[[10]]]
     diag(Sd) <- 1 - Sout$inbreeding

     if(!returnS) S <- NULL
     rm("Cout")

   } else{
stop("code not written to parallelize function") #FIXME
#        listA <- data.frame(Row = as.integer(rep(1:length(A@p[-1]), diff(A@p))), Column = as.integer(A@i + 1))
#        wrap_dij <- function(x){
#           sub_lA <- listA[min(x):max(x), 1:2]
#           lA_r <- dim(sub_lA)[1]
#           Cout <- .C("dijp",
#		as.integer(numeric.pedigree[, 2] - 1),
#		as.integer(numeric.pedigree[, 3] - 1), 
#                as.integer(lA_r),  
#		as.integer(sub_lA[, 1] - 1),  
#		as.integer(sub_lA[, 2] - 1),  
#		as.integer(A@i),  
#		as.integer(A@p), 
#		as.double(A@x/2),  
#		as.double(rep(0, lA_r)))
#         Cout[[9]]
#        }

#        cat("starting to make D...")
#        Dijs <- parallel::pvec(seq(1, dim(listA)[1], 1), FUN = wrap_dij, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE)
  
#        D <- Matrix(0, N, N, sparse = TRUE, dimnames = list(as.character(pedigree[, 1]), NULL))
#        D@uplo <- "U"
#        D@i <- A@i
#        D@p <- A@p
#        if(!returnA) A <- NULL
#        D@x <- Dijs
#        D <- drop0(D)
#        diag(D) <- 2 - dA

     }

  cat(".done", "\n")
  
  if(det) logDet <- determinant(Sd, logarithm = TRUE)$modulus[1] else logDet <- NULL
  if(invertSd){
    cat("starting to invert Sd...")
    Sdinv <- as(solve(Sd), "dgCMatrix")
      Sdinv@Dimnames <- Sd@Dimnames
    cat(".done", "\n")
    listSdinv <- sm2list(Sdinv, rownames=pedigree[,1], colnames=c("row", "column", "Dinverse"))
 return(list(S = S, Sd = Sd, logDet = logDet, Sdinv=Sdinv, listSdinv=listSdinv))
  } else{
    return(list(S = S, Sd = Sd, logDet = logDet))
    } 

#What is diagonal for heterogameitc sex?
## Should the matrix just be for those of 1 sex? Or can I automatically prune to 1 sex
### Model should probably just be for the 1 sex
## What about all forms of dosage compensation?


}


################################################################################
makeSd(FG90, heterogametic = "0", returnS = TRUE)

pedFS <- simPedHS(1, 1, 4)
makeSd(pedFS, heterogametic = "M", returnS = TRUE)





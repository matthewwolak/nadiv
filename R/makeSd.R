rm(list = ls())
library(nadiv)


makeSdR <- function(pedigree, heterogametic,
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
  N2 <- N + 1 #TODO necessary?
#FIXME turned off next check so can test parallel=TRUE on small pedigrees
#  if(parallel){
#    if(length(Sout$S@x)/ncores < 10){
#      warning("pedigree too small - 'parallel' set to FALSE instead")
#      parallel <- FALSE
#    }
#  }

  if(!parallel){
     cat("starting to make Sd...")

     Cout <- .C("sdij",
                as.integer(numeric.pedigree[, 2] - 1), 
		as.integer(numeric.pedigree[, 3] - 1), 
		as.integer(A@i), 			
		as.integer(A@p),                        
		as.double(A@x/2),                       
		as.integer(N),                           
		as.double(rep(0, length(A@x))),         
		as.integer(rep(0, length(A@i))),        
                as.integer(rep(0, N)),                  
		as.integer(0))	                        

     D <- Matrix(0, N, N, sparse = TRUE, dimnames = list(as.character(pedigree[, 1]), NULL))
     D@uplo <- "U"
     D@i <- Cout[[8]][1:Cout[[10]]]
     D@p <- c(Cout[[9]], Cout[[10]])
     D@x <- Cout[[7]][1:Cout[[10]]]
     diag(D) <- 2 - dA

     if(!returnA) A <- NULL
     rm("Cout")

   } else{
        listA <- data.frame(Row = as.integer(rep(1:length(A@p[-1]), diff(A@p))), Column = as.integer(A@i + 1))
        wrap_dij <- function(x){
           sub_lA <- listA[min(x):max(x), 1:2]
           lA_r <- dim(sub_lA)[1]
           Cout <- .C("dijp",
		as.integer(numeric.pedigree[, 2] - 1),
		as.integer(numeric.pedigree[, 3] - 1), 
                as.integer(lA_r),  
		as.integer(sub_lA[, 1] - 1),  
		as.integer(sub_lA[, 2] - 1),  
		as.integer(A@i),  
		as.integer(A@p), 
		as.double(A@x/2),  
		as.double(rep(0, lA_r)))
         Cout[[9]]
        }

        cat("starting to make D...")
        Dijs <- parallel::pvec(seq(1, dim(listA)[1], 1), FUN = wrap_dij, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE)
  
        D <- Matrix(0, N, N, sparse = TRUE, dimnames = list(as.character(pedigree[, 1]), NULL))
        D@uplo <- "U"
        D@i <- A@i
        D@p <- A@p
        if(!returnA) A <- NULL
        D@x <- Dijs
        D <- drop0(D)
        diag(D) <- 2 - dA

     }

  cat(".done", "\n")
  
  if(det) logDet <- determinant(D, logarithm = TRUE)$modulus[1] else logDet <- NULL
  if(invertD){
    cat("starting to invert D...")
    Dinv <- as(solve(D), "dgCMatrix")
      Dinv@Dimnames <- D@Dimnames
    cat(".done", "\n")
    listDinv <- sm2list(Dinv, rownames=pedigree[,1], colnames=c("row", "column", "Dinverse"))
 return(list(A = A, D = D, logDet = logDet, Dinv=Dinv, listDinv=listDinv))
  } else{
    return(list(A = A, D = D, logDet = logDet))
    } 

#What is diagonal for heterogameitc sex?
## Should the matrix just be for those of 1 sex? Or can I automatically prune to 1 sex
## What about all forms of dosage compensation?


}


################################################################################
makeSdR(FG90, heterogametic = "0", returnS = TRUE)




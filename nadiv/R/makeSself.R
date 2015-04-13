makeSself <- function(pedigree, heterogametic = "0", DosageComp = c(NULL, "ngdc", "hori", "hedo", "hoha", "hopi"), returnS = FALSE){

    if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")
    numped <- numPed(pedigree[, 1:3])

    damsex <- pedigree[unique(numped[, 2])[-1], 4]
    if(any(damsex == heterogametic)){
       pedname <- names(pedigree)
       pedigree <- pedigree[, c(1,3,2,4)]
       names(pedigree) <- pedname
       numped <- numPed(pedigree[, 1:3])
      warning("Assuming female heterogametic (e.g., ZZ/ZW) sex chromosome system")
    }
    sex <- rep(-998, dim(pedigree)[1])
    sex[homs <- which(pedigree[,4] != heterogametic)] <- 1
    sex[hets <- which(pedigree[,4] == heterogametic)] <- 0
    N <- dim(numped)[1]
    N2 <- N + 1

    dc.model <- match.arg(DosageComp)
    if(is.null(dc.model)){
      warning("Assuming 'ngdc' dosage compensation model")
      dc.model <- "ngdc"
    }


    if(dc.model == "ngdc"){
          dnmiss <- which(numped[,2] != -998)
          fsnmiss <- which(numped[,3] != -998 & sex == 1)
          bnmiss <- which(numped[, 2] != -998 & numped[, 3] != -998)
          nA <- N + 2 * length(dnmiss) + 2 * length(fsnmiss)
          nA <- nA + 2 * sum(duplicated(paste(numped[, 2], numped[, 3])[bnmiss]) == FALSE)
          Q.col <- c(numped[,1][dnmiss], numped[,1][fsnmiss], 1:N) 
          Q.row <- c(numped[,2][dnmiss], numped[,3][fsnmiss], 1:N)
          Q.x <- c(rep(-0.5, length(dnmiss)), rep(-1, length(fsnmiss)), rep(1, N))
          ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
          Q <- Matrix(0, N, N, sparse = TRUE)
          Q[1, 2] <- 1
          Q@i <- as.integer(Q.row[ord] -1)
          Q@p <- as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1)
          Q@x <- as.double(Q.x[ord])

          numped[numped == -998] <- N2
          Vii <- (sex + 1)/2
          f <- c(rep(0, N), -1)

          Cout <- .C("sinv",
	    as.integer(numped[, 2] - 1), #dam
	    as.integer(numped[, 3] - 1),  #sire
	    as.double(f),  #f
	    as.double(Vii),  #vii
            as.integer(Q@i),  #iQP
	    as.integer(c(Q@p, length(Q@x))),  #pQP
            as.double(Q@x),  #xQP
            as.integer(N),   #nQP
	    as.integer(length(Q@x)), #nzmaxQP	
	    as.integer(rep(0, nA)), #iSP
	    as.integer(rep(0, N + 1)), #pSP
	    as.double(rep(0, nA)), #xSP
	    as.integer(nA), #nzmaxSP
            as.double(0.25), #DC
            as.double(Vii))  #sex

          Sinv <- Matrix(0, N, N)
          Sinv[1, 2] <- 1
          Sinv@i <- Cout[[10]][1:Cout[[13]]]
          Sinv@p <- Cout[[11]]
          Sinv@x <- Cout[[12]][1:Cout[[13]]]
          Vii <- Cout[[4]]
          f <- Cout[[3]][1:N]



    } else{
         if(dc.model != "hopi"){
             fdnmiss <- which(numped[,2] != -998 & sex == 1)
             mdnmiss <- which(numped[,2] != -998 & sex == 0)
             fsnmiss <- which(numped[,3] != -998 & sex == 1)
             bnmiss <- which(numped[, 2] != -998 & numped[, 3] != -998)
             nA <- N + 2 * (length(fdnmiss) + length(mdnmiss)) + 2 * length(fsnmiss)
             nA <- nA + 2 * sum(duplicated(paste(numped[, 2], numped[, 3])[bnmiss]) == FALSE)
             Q.col <- c(numped[,1][fdnmiss], numped[,1][mdnmiss], numped[,1][fsnmiss], 1:N) 
             Q.row <- c(numped[,2][fdnmiss], numped[,2][mdnmiss], numped[,3][fsnmiss], 1:N)
             Q.x <- c(rep(-0.5, length(fdnmiss)), rep(-1, length(mdnmiss)), rep(-0.5, length(fsnmiss)), rep(1, N))
             ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
             Q <- Matrix(0, N, N, sparse = TRUE)
             Q[1, 2] <- 1
             Q@i <- as.integer(Q.row[ord] - 1)
             Q@p <- as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1)
             Q@x <- as.double(Q.x[ord])

             numped[numped == -998] <- N2
             Vii <- (2 - sex)
             f <- c(rep(0, N), -1)

             Cout <- .C("sinv",
	       as.integer(numped[, 2] - 1), #dam
	       as.integer(numped[, 3] - 1),  #sire
	       as.double(f),  #f
	       as.double(Vii),  #vii
               as.integer(Q@i),  #iQP
	       as.integer(c(Q@p, length(Q@x))),  #pQP
               as.double(Q@x),  #xQP
               as.integer(N),   #nQP
	       as.integer(length(Q@x)), #nzmaxQP	
	       as.integer(rep(0, nA)), #iSP
	       as.integer(rep(0, N + 1)), #pSP
	       as.double(rep(0, nA)), #xSP
	       as.integer(nA), #nzmaxSP
               as.double(1), #DC
               as.double(Vii))  #sex

             Sinv <- Matrix(0, N, N)
             Sinv[1, 2] <- 1
             Sinv@i <- Cout[[10]][1:Cout[[13]]]
             Sinv@p <- Cout[[11]]
             Sinv@x <- Cout[[12]][1:Cout[[13]]]
             Vii <- Cout[[4]]
             f <- Cout[[3]][1:N]


         } else{
                dnmiss <- which(numped[,2] != -998)
                nA <- N + 2 * length(dnmiss)
                nA <- nA + 2 * sum(duplicated(paste(numped[, 2], numped[, 3])[dnmiss]) == FALSE)
                Q.col <- c(numped[,1][dnmiss], 1:N) 
                Q.row <- c(numped[,2][dnmiss], 1:N)
                Q.x <- c(rep(-0.5, length(dnmiss)), rep(1, N))
                ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
                Q <- Matrix(0, N, N, sparse = TRUE)
                Q[1, 2] <- 1
                Q@i <- as.integer(Q.row[ord] - 1)
                Q@p <- as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1)
                Q@x <- as.double(Q.x[ord])

                numped[numped == -998] <- N2
                Vii <- rep(1, N) 
                f <- rep(0, N)

             Cout <- .C("sinv",
	       as.integer(numped[, 2] - 1), #dam
	       as.integer(numped[, 3] - 1),  #sire
	       as.double(f),  #f
	       as.double(Vii),  #vii
               as.integer(Q@i),  #iQP
	       as.integer(c(Q@p, length(Q@x))),  #pQP
               as.double(Q@x),  #xQP
               as.integer(N),   #nQP
	       as.integer(length(Q@x)), #nzmaxQP	
	       as.integer(rep(0, nA)), #iSP
	       as.integer(rep(0, N + 1)), #pSP
	       as.double(rep(0, nA)), #xSP
	       as.integer(nA), #nzmaxSP
               as.double(0), #DC
               as.double(Vii))  #sex

             Sinv <- Matrix(0, N, N)
             Sinv[1, 2] <- 1
             Sinv@i <- Cout[[10]][1:Cout[[13]]]
             Sinv@p <- Cout[[11]]
             Sinv@x <- Cout[[12]][1:Cout[[13]]]
             Vii <- Cout[[4]]
             f <- Cout[[3]][1:N]

           }
      }              


    listSinv <- sm2list(Sinv, rownames = as.character(pedigree[, 1]), colnames = c("Row", "Column", "Sinverse"))
    Sinv@Dimnames <- list(pedigree[,1], NULL) 

    if(returnS){
       cat(paste("S-inverse made: Starting to make S..."))
          T <- as(solve(Q), "dgCMatrix")
          S <- as(t(T) %*% Diagonal(N, Vii) %*% T, "dgCMatrix")
       cat(paste(".done", "\n"))
    } else{
         S <- NULL
      }

return(list(model = dc.model, S = S, Sinv = Sinv, listSinv = listSinv, inbreeding = f, v = Vii))
}


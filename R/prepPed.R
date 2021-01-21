#' Prepares a pedigree by sorting and adding 'founders'
#' 
#' This function takes a pedigree, adds missing founders, and then sorts the
#' pedigree.
#' 
#' Many functions (both in nadiv and from other programs) dealing with
#' pedigrees must first sort a pedigree such that individuals appear in the ID
#' column in rows preceding where they appear in either the Dam or Sire
#' column.  Further, these functions and programs require that all individuals
#' in the dam and sire columns of a pedigree also have an entry in the ID
#' column.  This function easily prepares data sets to accommodate these
#' requirements using a very fast topological sorting algorithm.
#' 
#' NOTE: more columns than just a pedigree can be passed in the \code{pedigree}
#' argument.  In the case of missing founders, these columns are given NA
#' values for all rows where founders have been added to the pedigree.  The
#' entire object supplied to \code{pedigree} is ordered, ensuring that all
#' information remains connected to the individual
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' When a non-null argument is given to \code{gender}, dams without an entry in
#' the ID column (that are subsequently added to the pedigree) are given the
#' gender designated for other dams (and similarly for sires).
#' 
#' The \code{check} argument performs checks on the format of the pedigree
#' supplied to try and identify any issues regarding the notation of missing
#' values and validity of the basic pedigree for further processing.
#' 
#' @param pedigree An object, where the first 3 columns correspond to: ID, Dam,
#'   & Sire. See details.
#' @param gender An optional character for the name of the column in
#'   \code{pedigree} that corresponds to the gender/sex of individuals. If
#'   specified, \code{prepPed} will assign a gender to any founders it adds to
#'   the pedigree.
#' @param check A logical argument indicating if checks on the validity of the
#'   pedigree structure should be made
#'
#' @return The pedigree object (can have more columns than just ID, Dam, and
#' Sire), where: (1) the ID column contains an ID for all individuals from the
#' original pedigree object's ID, Dam, and Sire columns (i.e., founders are
#' added) and (2) the pedigree is now sorted so that individuals are not in
#' rows preceding either their Dam or Sire.
#' @seealso \code{\link[nadiv]{genAssign}}, \code{\link[nadiv]{prunePed}}
#' @examples
#' 
#' # First create an unordered pedigree with (4) missing founders
#'   warcolak_unsuitable <- warcolak[sample(seq(5, nrow(warcolak), 1),
#' 	size = (nrow(warcolak) - 4), replace = FALSE), ]
#'   nrow(warcolak)
#'   nrow(warcolak_unsuitable)
#' # Fix and sort the pedigree
#' ## Automatically assign the correct gender to the added founders
#' ### Also sort the data accompanying each individual
#'   warcolak_fixed_ordered <- prepPed(warcolak_unsuitable, gender = "sex")
#'   head(warcolak_fixed_ordered)
#' 
#' @export
prepPed <- function(pedigree, gender = NULL, check = TRUE){

 self <- FALSE
 if(check){      
   if(length(d0 <- which(pedigree[, 2] == 0)) > 0){
     pedigree[d0, 2] <- NA
     warning("Zero in the dam column interpreted as a missing parent")
   }
   if(length(s0 <- which(pedigree[, 3] == 0)) > 0){
     pedigree[s0, 3] <- NA
     warning("Zero in the sire column interpreted as a missing parent")
   }
   if(length(dast <- which(pedigree[,2] == "*")) > 0) pedigree[dast, 2] <- NA
   if(length(sast <- which(pedigree[,3] == "*")) > 0) pedigree[sast, 3] <- NA
   if(all(is.na(pedigree[, 2])) & all(is.na(pedigree[, 3]))){
     stop("All dams and sires are missing")
   }
   if(any(idMiss <- pedigree[, 1] == 0 | pedigree[, 1] == "0" | pedigree[, 1] == "*" | is.na(pedigree[, 1]))){
     warning("Missing value in the ID column - row discarded")
     warning("Check to ensure first three columns of the pedigree object are ID, Dam, and Sire")
     pedigree <- pedigree[-which(idMiss), ]
   }

   if(any(pedigree[!is.na(pedigree[, 2]), 2] %in% pedigree[!is.na(pedigree[, 3]), 3])){
    self <- TRUE
    warning("Dams appearing as Sires - assumed selfing in pedigree")
   }

 }  #<-- end if check

  
 facflag <- is.factor(pedigree[, 1])
 udam <- if(facflag) as.character(unique(pedigree[, 2])) else unique(pedigree[, 2])
 udam <- udam[!is.na(udam)]
 missdam <- udam[which(is.na(match(udam, pedigree[, 1])))]
 usire <- if(facflag) as.character(unique(pedigree[, 3])) else unique(pedigree[, 3])
 usire <- usire[!is.na(usire)]
 misssire <- usire[which(is.na(match(usire, pedigree[, 1])))]


 if(length(missdam) == 0 & length(misssire) == 0){
   ped_fixed <- pedigree
 } else{
     missparent <- union(missdam, misssire)
     topPed <- data.frame(missparent,
         rep(NA, length(missparent)),
         rep(NA, length(missparent)),
         matrix(NA, nrow = length(missparent), ncol = ncol(pedigree) - 3))
     names(topPed) <- names(pedigree)
     if(!is.null(gender)){
        if(self | (length(missdam) + length(misssire)) != length(missparent)){
          warning("Selfing in the pedigree: the gender argument might not apply")
        }
        if(is.factor(pedigree[, gender])){
          damgender <- as.character(pedigree[which(pedigree[, 1] == udam[which(!udam %in% missdam)][1]), gender]) 
          siregender <- as.character(pedigree[which(pedigree[, 1] == usire[which(!usire %in% misssire)][1]), gender]) 
        } else{
             damgender <- pedigree[which(pedigree[, 1] == udam[which(!udam %in% missdam)][1]), gender] 
             siregender <- pedigree[which(pedigree[, 1] == usire[which(!usire %in% misssire)][1]), gender]
          } 
        topPed[, gender] <- c(rep(damgender, length(missdam)), rep(siregender, length(misssire)))
     }
     ped_fixed <- rbind(topPed, pedigree)
   }
 npf <- nrow(ped_fixed)

 if(sum((na.omit(ped_fixed[, 2]) %in% ped_fixed[, 1]) == FALSE) > 0 & any(is.na(ped_fixed[, 2]) == FALSE)){
   stop("Something wicked happened: individuals appearing as dams but not added to pedigree. Report as possible bug to <matthewwolak@gmail.com>")
 }
 if(sum((na.omit(ped_fixed[, 3]) %in% ped_fixed[, 1]) == FALSE) > 0 & any(is.na(ped_fixed[, 3]) == FALSE)){
   stop("Something wicked happened: individuals appearing as sires but not added to pedigree. Report as possible bug to <matthewwolak@gmail.com>")
 }
 if(sum(duplicated(ped_fixed[, 1])) > 0){
   stop("some individuals appear more than once in the pedigree")
 }

 nPed_fixed <- numPed(ped_fixed[, 1:3], check = FALSE)
 Cout <- .C("gaUnsort", PACKAGE = "nadiv",
	as.integer(nPed_fixed[, 2] - 1),
	as.integer(nPed_fixed[, 3] - 1),
        as.integer(rep(0, npf)),
	as.integer(rep(0, npf)),
	as.integer(npf))

 # check for errors/return values indicating an infinite loop
 if(Cout[[5]] < npf){
   infLooper <- Cout[[5]] + 1
   if(Cout[[3]][1] == -999){
     cat("\n Problem: infinite pedigree loop involving individual",
       as.character(ped_fixed[infLooper, 1]), "and dam/female ancestors\n\n")
   }
   if(Cout[[4]][1] == -999){
     cat("\n Problem: infinite pedigree loop involving individual",
       as.character(ped_fixed[infLooper, 1]), "and sire/male ancestors\n\n")
   }
   stop("Check for individuals that appear as their own ancestor")
 }

  # Order by number of dam/sire paths to dam/sire founder:
  ## first by parent paths wit MAX then by parent with MIN
  ped_fixed_ord <- ped_fixed[order(pmax.int(Cout[[3]], Cout[[4]]),
  				    pmin.int(Cout[[3]], Cout[[4]])), ]
  				    
  # Make this into newly ordered numPed (replace existing object)
  nPed_fixed <- numPed(ped_fixed_ord, check = FALSE)
  generation <- rep(-1, npf)  
    generation[which(nPed_fixed[, 2] == -998 & nPed_fixed[, 3] == -998)] <- 0
   
  Cout <- .C("ga", PACKAGE = "nadiv",
	as.integer(nPed_fixed[, 2] - 1),
	as.integer(nPed_fixed[, 3] - 1),
        as.integer(generation),
	as.integer(npf))   
  generation[] <- Cout[[3]]    				    
  ped_fixed_ord <- ped_fixed_ord[order(generation), ]
  itwork <- try(expr = numPed(ped_fixed_ord[, 1:3]))

 return(ped_fixed_ord) 
}


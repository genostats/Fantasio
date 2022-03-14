#' Computation of FLOD scores
#' 
#' This function is uses to compute FLOD scores on individual present in a sample
#' 
#' @param submaps a atlas object
#' @param w.id The indices of individals for which FLOD are computed
#' @param q The assumed frequency of the mutation involved in the disease for each individual. Default is 0.0001.
#' 
#' @details This function iterates over the slots submaps_list of the atlas object.
#' @details For each submaps in the slots submaps_list of the object, the slot FLOD.prob will be filled with a matrix of dimension : number_inidividual x number_of_markers
#' 
#' @return the atlas object with each FLOD.prob slot of each submaps in the slot submaps_list computed
#' 
#' @seealso setHBDProb
#' @seealso setHFLOD
#' 
#' 
#' @export
setFLOD <- function(submaps, w.id, q = 1e-4) {
  if(class(submaps@submaps_list[[1]])[1] != "snpsMatrix" & class(submaps@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  moy_FLOD <- matrix(0, nrow = length(w.id), ncol = atlas@submaps_list[[1]]@ncol)#HBD matrix
  dimnames(moy_FLOD) <- list( rownames(submaps@submaps_list[[i]]@HBD.prob) , lapply(c(1:atlas@submaps_list[[1]]@ncol), function(i) paste0('s', i)) )
  
  #Computation of FLOD score with the formula
  for(i in seq_along(submaps@submaps_list))
  {
    FLOD <- matrix(0.0, nrow = nrow(submaps@submaps_list[[i]]@HBD.prob), ncol = submaps@submaps_list[[i]]@ncol)
    rownames(FLOD) <- rownames(submaps@submaps_list[[i]]@HBD.prob) 
    colnames(FLOD) <- colnames(submaps@submaps_list[[i]]@HBD.prob)
    
    for (j in seq_len(nrow(FLOD))) {
      if( (submaps@submaps_list[[i]]@a[w.id[j]] < 1) & !is.na(submaps@submaps_list[[i]]@f[w.id[j]]) ) {
        FLOD[j,] <- log10((submaps@submaps_list[[i]]@HBD.prob[j,] + q * ( 1 - submaps@submaps_list[[i]]@HBD.prob[j,]))/
                                 (submaps@submaps_list[[i]]@f[w.id[j]] + q * ( 1 - submaps@submaps_list[[i]]@f[w.id[j]]))) 
      }
    }
    moy_FLOD <- moy_FLOD + FLOD/length(submaps@submaps_list)
  }
  submaps@FLOD_recap <- moy_FLOD
  return(submaps)
}


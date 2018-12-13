#' Computation of FLOD scores
#' 
#' This function is uses to compute FLOD scores on individual present in a sample
#' 
#' @param submaps a atlas object
#' @param condition A vector containing a list of individuals. Only this list of individuals will have their FLOD scores computed
#' @param q Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual. Default is 0.0001.
#' 
#' @details This function iterates over the slots submaps_list of the atlas object.
#' @details For each submaps in the slots submaps_list of the object, the slot FLOD.prob will be filled with a matrix of dimension : number_inidividual x number_of_markers
#' 
#' @return the atlas object with each FLOD.prob slot of each submaps in the slot submaps_list computed
#' 
#' @seealso setHBDprob
#' @seealso setHFLOD
#' 
#' 
#' @export
setFLOD <- function(submaps, condition, q = 1e-4)
{
  if(class(submaps@submaps_list[[1]])[1] != "snpsMatrix" & class(submaps@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  #Computation of FLOD score with the formula
  for(i in 1:length(submaps@submaps_list))
  {
    FLOD_prob <- matrix(0.0, nrow = nrow(submaps@submaps_list[[i]]@HBD.prob), ncol = submaps@submaps_list[[i]]@ncol)
    dimnames(FLOD_prob) <- list(rownames(FLOD_prob) <- rownames(submaps@submaps_list[[i]]@HBD.prob), colnames(FLOD_prob) <- colnames(submaps@submaps_list[[i]]@HBD.prob))
    
    for (j in 1:nrow(FLOD_prob))
    {
      if( (submaps@submaps_list[[i]]@a[condition[j]] < 1) & !is.na(submaps@submaps_list[[i]]@f[condition[j]]) )
      {
        FLOD_prob[j,] <- log10((submaps@submaps_list[[i]]@HBD.prob[j,] + q * ( 1 - submaps@submaps_list[[i]]@HBD.prob[j,]))/
                                 (submaps@submaps_list[[i]]@f[condition[j]] + q * ( 1 - submaps@submaps_list[[i]]@f[condition[j]]))) 
      }
    }
    submaps@submaps_list[[i]]@FLOD <- FLOD_prob 
  }
  return(submaps)
  
}

#' Computation of FLOD scores
#' 
#' This function is uses to compute FLOD scores on individual present in a sample
#' 
#' @param submaps a list.submaps object
#' @param condition A vector containing a list of individuals. Only this list of individuals will have their FLOD scores computed
#' @param q Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual. Default is 0.0001.
#' 
#' @details This function iterates over the slots atlas of the list.submaps object.
#' @details For each submaps in the slots atlas of the object, the slot FLOD.prob will be filled with a matrix of dimension : number_inidividual x number_of_markers
#' 
#' @return the list.submaps object with each FLOD.prob slot of each submaps in the slot atlas computed
#' 
#' @seealso set.HBD.prob
#' @seealso set.HFLOD
#' 
#' 
#' @export
set.FLOD <- function(submaps, condition, q = 1e-4)
{
  if(class(submaps@atlas[[1]])[1] != "snps.matrix" & class(submaps@atlas[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps.") 
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  #Computation of FLOD score with the formula
  for(i in 1:length(submaps@atlas))
  {
    FLOD_prob <- matrix(0.0, nrow = nrow(submaps@atlas[[i]]@HBD.prob), ncol = submaps@atlas[[i]]@ncol)
    dimnames(FLOD_prob) <- list(rownames(FLOD_prob) <- rownames(submaps@atlas[[i]]@HBD.prob), colnames(FLOD_prob) <- colnames(submaps@atlas[[i]]@HBD.prob))
    
    for (j in 1:nrow(FLOD_prob))
    {
      if( (submaps@atlas[[i]]@a[condition[j]] < 1) & !is.na(submaps@atlas[[i]]@f[condition[j]]) )
      {
        FLOD_prob[j,] <- log10((submaps@atlas[[i]]@HBD.prob[j,] + q * ( 1 - submaps@atlas[[i]]@HBD.prob[j,]))/
                                 (submaps@atlas[[i]]@f[condition[j]] + q * ( 1 - submaps@atlas[[i]]@f[condition[j]]))) 
      }
    }
    submaps@atlas[[i]]@FLOD <- FLOD_prob 
  }
  return(submaps)
  
}

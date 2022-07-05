#' HBD/FLOD recap 
#' 
#' This function creates HBD and FLOD recap dataframe.
#' 
#' @param atlas a atlas object
#' @param recap.by.segments : whether you want the recap by segments or snps
#' 
#' @details For each individual and each marker the function computes 
#' @details the mean value of every HBD probabilities computed through the submaps.
#' @details The same is done with the FLOD scores.  
#' @details This function returns a list of two dataframes with HBD and FLOD in it.
#' 
#' @seealso Fantasio
#' @seealso makeAtlasByDistance
#' @seealso segmentsListByHotspots
#' @seealso festim
#' @seealso setHBDProb
#' @seealso setFLOD
#' @seealso setHFLOD
#' @seealso HBDSegments
#'
#' @return This function returns a list of dataframes. 
#'
#' @export
recap <- function(atlas, recap.by.segments = atlas@bySegments) {
  
  if(class(atlas@submaps_list[[1]])[1] != "snpsMatrix" & class(atlas@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  if(class(atlas@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
    
  if(length(atlas@submaps_list) == 1) {
    matrice_HBD <- atlas@submaps_list[[1]]@HBD.prob
    matrice_FLOD <- atlas@submaps_list[[1]]@FLOD
    l <- list(matrice_HBD, matrice_FLOD)
    return(l)	
  }
  
  # HBD probabilities and FLOD scores
  proba_HBD  <- submapHBDRecap(atlas)
  FLOD <- submapFLODRecap(atlas)
  
  if(recap.by.segments) {
    recapBySegments(atlas, proba_HBD, FLOD)
  } else {
    recapBySnps(atlas, proba_HBD, FLOD)
  }
}  

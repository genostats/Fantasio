#' HBD/FLOD recap 
#' 
#' This function creates HBD and FLOD recap dataframe.
#' 
#' @param submaps a atlas object
#' @param recap.by.segments : whether you want the recap by segments or snps
#' @param list.id : a list individual 
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
#' @seealso setHBDprob
#' @seealso setFLOD
#' @seealso setHFLOD
#' @seealso HBDsegments
#'
#' @return This function returns a list of dataframes. 
#'
#' @export
recap <- function(submaps, recap.by.segments = submaps@bySegments, list.id) {
  
  if(class(submaps@submaps_list[[1]])[1] != "snpsMatrix" & class(submaps@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
    
  if(!missing(list.id))
  {
    if(list.id == "all")
    {
      condition <- 1:nrow(submaps@submap_summary)
    }else{
      vec <- strsplit(list.id, "_")
      condition <- sapply(vec, function(i) which(submaps@submap_summary$FID == i[1] & submaps@submap_summary$IID == i[2]))
    }
  }else{
    condition <- which(submaps@submap_summary$QUALITY >= 95 & submaps@submap_summary$INBRED)
  }
  
  if(length(submaps@submaps_list) == 1)
  {
  	matrice_HBD <- submaps@submaps_list[[1]]@HBD.prob
  	matrice_FLOD <- submaps@submaps_list[[1]]@FLOD
  	l <- list(matrice_HBD, matrice_FLOD)
  	return(l)	
  }
  
  #list form for HBD probabilities and FLOD scores
  proba_HBD  <- submap.HBD(submaps)
  proba_FLOD <- submap.FLOD(submaps)
  
  if(recap.by.segments)
  {
    recap.by.segments(submaps, proba_HBD, proba_FLOD)
  }
  else
  {
    recap.by.snps(submaps, proba_HBD, proba_FLOD)
  }
}  

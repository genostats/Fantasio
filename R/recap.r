#' HBD/FLOD recap 
#' 
#' This function creates HBD and FLOD recap dataframe.
#' 
#' @param submaps a list.submaps object
#' @param recap.by.segments : whether you want the recap by segments or snps
#' @param list.id : a list individual 
#' 
#' @details For each individual and each marker the function computes 
#' @details the mean value of every HBD probabilities computed through the submaps.
#' @details The same is done with the FLOD scores.  
#' @details This function returns a list of two dataframes with HBD and FLOD in it.
#' 
#' @seealso Fantasio
#' @seealso makeSubmapsBySnps
#' @seealso createSegmentsListByHotspots
#' @seealso festim
#' @seealso set.HBD.prob
#' @seealso set.FLOD
#' @seealso set.HFLOD
#' @seealso HBD.segments
#'
#' @return This function returns a list of dataframes. 
#'
#' @export
recap <- function(submaps, recap.by.segments=FALSE, list.id)
{
  
  if(class(submaps@atlas[[1]])[1] != "snps.matrix" & class(submaps@atlas[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.") 
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix to eat")
    
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
  
  if(length(submaps@atlas) == 1)
  {
  	matrice_HBD <- submaps@atlas[[1]]@HBD.prob
  	matrice_FLOD <- submaps@atlas[[1]]@FLOD
  	l <- list(matrice_HBD, matrice_FLOD)
  	return(l)	
  }
  
  #recuperer les probas HBD pour chaque sous cartes sous forme d'une liste
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

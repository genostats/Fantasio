#' Number of marker selected in a submap
#' 
#' This function gives for each marker the number of time that it has been selected in a submap. 
#' 
#' @param submaps a list of submaps 
#' @param bedmatrix a bed.matrix object
#' 
#' @details the dataframe contains the following elements :
#' @details - by columns is the number of markers used
#' @details - by row is the number of time picked
#' 
#' @return this function returns a dataframe.
#' 
#' @seealso setSummary
#' 
#' @export
summaryMarker <- function(submaps, bedmatrix)
{
  if(class(submaps[[1]])[1] != "snps.matrix" & class(submaps[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.") 
  if(class(bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix to eat")
  
  
  b <- summaryMap(submaps)
  res <- c()
  for(i in 1:length(submaps))
  {
    taille <- length(which(b$Freq == i))
    if(taille == 0)
      break()
    
    res <- c(res, taille)
      
    
  }
  zero <- length(bedmatrix@snps$chr) - sum(res)
  
  df <- data.frame(
    number_of_time_picked = 0:length(res),
    number_of_markers = c(zero, res)
  )
  df
}
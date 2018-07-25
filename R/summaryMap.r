#' Number of selected marker through the submaps
#' 
#' This function  will give the number of markers selected in your submaps.
#' 
#' @param submaps a list containing the different submaps created
#' 
#' @details make sure to provide a list of submaps in argument. This function returns a dataframe with 2 columns : 
#' @details -snsp : the id of the snps 
#' @details -Freq : the number of times the marke has been picked
#' 
#' @return  This function returns a dataframe.
#' 
#' @seealso setSummary
#' 
#' @export
summaryMap <- function(submaps)
{
  if(class(submaps[[1]])[1] != "snps.matrix" & class(submaps[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps.") 
  
  snps <- unlist(sapply(submaps, function(x) x@map$id)) #every markers that has been choosen on submaps
  b <- as.data.frame(table(snps),  stringsAsFactors=FALSE)
  b
}


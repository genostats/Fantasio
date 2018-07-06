#' Number of marker selected in a submap
#' 
#' This function gives for each submap the number of marker selected. 
#' 
#' @param submaps a list.submaps object
#' 
#' @details the dataframe contains the following elements :
#' @details - by columns is the number of markers used
#' @details - by row is each submap
#' 
#' @return this function returns a dataframe.
#' 
#' @seealso \code{\link{}}
#' 
#' @examples  
#' 
#' @export
markerSummary <- function(h)
{
  df <- as.data.frame(sapply(h@atlas, function(x) ncol(x)))
  n_submaps <- c()
  for(i in 1:length(h@atlas)) n_submaps <- c(n_submaps, paste("Submap", i))
  rownames(df) <- n_submaps
  colnames(df) <- "number_of_markers_used"
  df
}


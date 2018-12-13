#' Number of marker selected in a submap
#' 
#' This function gives for each submap the number of marker selected. 
#' 
#' @param submaps a atlas object
#' 
#' @details the dataframe contains the following elements :
#' @details - by columns is the number of markers used
#' @details - by row is each submap
#' 
#' @return this function returns a dataframe.
#' 
#' @seealso setSummary
#' 
#' 
#' @export
markerSummary <- function(submaps)
{
  df <- as.data.frame(sapply(submaps@submaps_list, function(x) ncol(x)))
  n_submaps <- numeric(length(submaps@submaps_list))
  for(i in 1:length(submaps@submaps_list)) n_submaps[i] <- paste("Submap", i)
  rownames(df) <- n_submaps
  colnames(df) <- "number_of_markers_used"
  df
}


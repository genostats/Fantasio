#'Find HBD segments
#' 
#'This function creates HBD segments used for creating plots                      
#'
#' @param submaps a list.submaps object
#' @param n.consecutive.marker the number of consecutive marker with a probabilitie equal or greater to the value of the threshold, to be use to fing HBD segments (default is 5)
#' @param threshold the minimum value of HBD probabilities given for a marker (default is 0.5)
#' @param recap.by.segments if the summary files has to be computed considering segments or snps (defaut is FALSE) 
#' 
#' @details The threshold is the minimum value from which we consider the marker is HBD. From this marker we want a minumum number of consecutive markers to create a HBD segment.
#' (argument : n.consecutive.marker)
#' 
#' @return This function returns a list of dataframe with 11 columns : 
#' @return - individual
#' @return - family
#' @return - status
#' @return - start
#' @return - end
#' @return - size
#' @return - chromosome
#' @return - start_pos
#' @return - end_pos
#' @return - start_dist
#' @return - end_dist
#' 
#' @seealso setSummary
#' @seealso set.HBD.prob
#' 
#' @examples  
#' 
#' @export
HBD.segments <- function(submaps, n.consecutive.marker = 5, threshold = 0.5, recap.by.segments=FALSE)
{

  if(recap.by.segments)
  {
    HBD.segments.by.segments(submaps, submaps@HBD_recap, n.consecutive.marker, threshold)
  }
  else
  {
    HBD.segments.by.snps(submaps, submaps@HBD_recap, n.consecutive.marker, threshold)
  }
}


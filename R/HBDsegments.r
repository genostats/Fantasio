#'Find HBDsegments
#' 
#'This function creates HBDsegments used for creating plots                      
#'
#' @param submaps a atlas object
#' @param n.consecutive.marker the number of consecutive marker with a probabilitie equal or greater to the value of the threshold, to be use to fing HBDsegments (default is 5)
#' @param threshold the minimum value of HBD probabilities given for a marker (default is 0.5)
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
#' @seealso setHBDprob
#' 
#' 
#' @export
HBDsegments <- function(submaps, n.consecutive.marker = 5, threshold = 0.5)
{
  
  if(submaps@bySegments)
  {
    HBDsegments.by.segments(submaps, submaps@HBD_recap, n.consecutive.marker, threshold)
  }
  else
  {
    HBDsegments.by.snps(submaps, submaps@HBD_recap, n.consecutive.marker, threshold)
  }
}


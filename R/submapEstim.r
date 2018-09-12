#'Summury of estimation through the submaps
#'
#'This function creates a summary of the a and f estimationq on submaps created. 
#' 
#' @param submaps a list of submaps
#' 
#' @details The first element of the returned list is the estimation of f through the submaps. 
#' @details The second element of the returned list is the estimation of a through the submaps. 
#' 
#' @return this fonction returns a list of dataframe.
#' 
#' @seealso setSummary
#' 
#' 
#' @export
submapEstim <- function(submaps)
{
  if(class(submaps[[1]])[1] != "snpsMatrix" & class(submaps[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  estimation <- list()
  dfA <- data.frame("estimation of a.submap" = sapply(submaps, function(x) x@a))
  dfF <- data.frame("estimation of f.submap" = sapply(submaps, function(x) x@f))
  
  estimation[[1]] <- dfF
  estimation[[2]] <- dfA
   
  return(estimation)
}
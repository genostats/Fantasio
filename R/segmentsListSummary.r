#' Summary of marker picked
#'
#' This function returns the number of segments and number of marker in the list of segments outputed by the function 
#' `segmentsListByHotspots`
#'
#' @param segmentList A list of segments
#'
#' @examples
#' #Please refer to vignette 
#'
#'
#' @export
segmentsListSummary <- function(segmentList)
{
  if(class(segmentList)[1] == "snpsSegments")
    segmentList <- segmentList@snpsSegments
  else if(class(segmentList)[1] != "HostspotsSegments" & class(segmentList)[1] != "list")
    stop("Argument must be of class 'snpsSegments', 'HostspotsSegments' or 'list'")
    
  #number of segments
  n_seg <- numeric(length(segmentList))

  for(i in seq_along(segmentList))
    n_seg[i] <- length(segmentList[[i]])
  #n_seg[23] <- NA
  #n_seg <- n_seg[!is.na(n_seg)]
  
  
  #number of markers 
  n_mark <- numeric(length(segmentList))
  for(i in seq_along(segmentList))
  {
    res <- numeric(length(segmentList[[i]]))
    for(j in seq_along(segmentList[[i]]))
    {
      res[j] <- length(segmentList[[i]][[j]])
    }
    res <- sum(res)
    n_mark[i] <- res
  }
  
  
  #dataframe
  df <- data.frame(
    chromosome = if(class(segmentList)[1] == "HostspotsSegments") seq_along(segmentList) else getOption("gaston.autosomes"),
    number_of_segments = n_seg, 
    number_of_markers= n_mark
  )
  df
}


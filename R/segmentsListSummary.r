#' Summary of marker picked
#'
#' This function is uses to return the number of segments and number of marker in the list of segments outputed by the function 
#' `createSegmentsListByHotspots`
#'
#' @export
segmentsListSummary <- function(segmentList)
{
  #number of segments
  n_seg <- c()
  for(i in 1:length(segmentList))
  {
    n_seg <- c(n_seg, length(segmentList[[i]]))
  }
  n_seg[23] <- NA
  n_seg <- n_seg[!is.na(n_seg)]
  
  
  #number of markers 
  n_mark <- c()
  for(i in 1:length(segmentList))
  {
    res <- c()
    for(j in 1:length(segmentList[[i]]))
    {
      res <- c(res, length(segmentList[[i]][[j]]))
    }
    res <- sum(res)
    n_mark <- c(n_mark, res)
  }
  
  
  #dataframe
  df <- data.frame(
    chromosom = 1:22,
    number_of_segments = n_seg, 
    number_of_markers= n_mark
  )
  df
}


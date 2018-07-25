#' Summary of marker picked
#'
#' This function is uses to return the number of segments and number of marker in the list of segments outputed by the function 
#' `createSegmentsListByHotspots`
#'
#' @param segmentList A list of segments
#'
#' @examples
#' ##install.packages("HGDP.CEPH", repos="https://genostats.github.io/R/") ## make this only one time
#' require(Fantasio)
#' require(HGDP.CEPH)
#' filepath <-system.file("extdata", "hgdp_ceph.bed", package="HGDP.CEPH")
#' x <- read.bed.matrix(filepath)
#' x <- set.stats(x)
#' x.me <- select.inds(x, population == "Bedouin")
#' x.me@ped$pheno <- rep(2,48) #The package analyzes only individualw with a status of 2
#' s <- createSegmentsListByHotspots(x.me)
#' segmentsListSummary(s)
#'
#' @export
segmentsListSummary <- function(segmentList)
{
  if(class(segmentList)[1] != "hotspot.segments")
    segmentList <- segmentList@snps.segments
    
  #number of segments
  n_seg <- numeric(length(segmentList))

  for(i in 1:length(segmentList))
  {
    n_seg[i] <- length(segmentList[[i]])
  }
  #n_seg[23] <- NA
  #n_seg <- n_seg[!is.na(n_seg)]
  
  
  #number of markers 
  n_mark <- numeric(length(segmentList))
  for(i in 1:length(segmentList))
  {
    res <- numeric(length(segmentList[[i]]))
    for(j in 1:length(segmentList[[i]]))
    {
      res[j] <- length(segmentList[[i]][[j]])
    }
    res <- sum(res)
    n_mark[i] <- res
  }
  
  
  #dataframe
  df <- data.frame(
    chromosome = if(class(segmentList)[1] == "hotspot.segments") 1:length(segmentList) else getOption("gaston.autosomes"),
    number_of_segments = n_seg, 
    number_of_markers= n_mark
  )
  df
}


#' Creation of a submaps
#' 
#' This function creates a submap using the list segments created by using gap between markers
#' 
#' @param bedmatrix a bed.matrix object 
#' @param segmentsList a list of segments for each chromosome
#' @param epsilon genotype error rate (default is 0.001)
#' @param snpIndices (optional) You have the possibility to pass your own list of markers
#' 
#' @details If snpIndices is given, the function creates a submap corresponding to the given SNPs.
#' @details Otherwise, the function iterates over the list of segments, then for each segments it picks randomly one marker 
#' and from this marker choses a marker every "step" step from back to end.
#' 
#' @return return an snpsMatrix object 
#' 
#' @seealso makeAllSubmapsByDistance
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' 
#' @export
createSubmapByDistance <- function(bedmatrix, segmentsList, epsilon = 1e-3, snpIndices)
{
  if(class(segmentsList)[1] != "snpsSegments")
    stop("mismatch segments list, need a list of segments created by the function 'segmentsListByDistance' ")
  
  unit <- segmentsList@unit
  step <- segmentsList@gap 
  
  if( unit != "Bases" & unit != "cM")
    stop("Error only cM or Bases are accepted, please make sure to choose between 'cM' or 'Bases' when creating segments list.")
  
    
  
  if(!missing(snpIndices))
  {
    submap <- snpIndices 
  } else {
    submap <- c()
    for(chr in 1:length(segmentsList@snpsSegments))
    {
      map <- segmentsList@snpsSegments[[chr]]
      if(length(map) > 0)
      {
        v <- getMarkerChromosomByDistance(x=bedmatrix, map=map, pas=step, unit=unit) #return an index vector of the marker pick randomly in the segment
        
        if(unit == "Bases")
        {
          tmp_dist <- bedmatrix@snps$pos[v]
        }
        else{
          tmp_dist <- bedmatrix@snps$dist[v]
        }
        #check the distance between markers in the end of two mini segments
        tmp <- diff(tmp_dist)
        v <- v[which(tmp >= step)]
        submap <- c(submap, v)
      }
    }
  }
  
  
  map <- bedmatrix@snps[submap , c("id","chr")]
  if(unit == "cM")
  {
    if(all(bedmatrix@snps$dist == 0)) {
      map$distance <- bedmatrix@snps$pos[submap]*1e-6
    } else {
      map$distance <- bedmatrix@snps$dist[submap]
    }
  }
  else
  {
    if(all(bedmatrix@snps$pos == 0)) {
      map$distance <- bedmatrix@snps$dist[submap]*1e6
    } else {
      map$distance <- bedmatrix@snps$pos[submap]
    }
  }
  
  
  log.emiss <- bed.logEmiss(bedmatrix, map=submap, epsilon=epsilon)
  
  new("snpsMatrix", step, length(submap), nrow(bedmatrix), submap, 
      bedmatrix@ped[,c("famid", "id", "father", "mother", "sex", "pheno")],
      map, log.emiss, epsilon)
}

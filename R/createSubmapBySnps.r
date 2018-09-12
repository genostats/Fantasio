#' Creation of a submaps
#' 
#' This function creates a submaps using the list segments created by using gap between markers
#' 
#' @param bedmatrix a bed.matrix object 
#' @param segmentsList a list of segment for each chromosomes
#' @param epsilon genotype error rate (default is 0.001)
#' @param fileName You have the possibility to pass your own list of marker, that is to say, your own submaps 
#' 
#' @details This function will iterates over the list of segments, then for each segments it will pick randomly one marker 
#' and from this marker chose a marker every "step" step from back to end.
#' @details Once the iteration over the list is over, the function will create an object and filled some of his slot.
#' @details If you are using a fileName, please make sure to have one marker per line
#' 
#' @return return an snpsMatrix object with some of his slots filled.
#' 
#' @seealso makeSubmapsBySnps
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' 
#' @export
createSubmapBySnps <- function(bedmatrix, segmentsList, epsilon = 1e-3, fileName)
{
  if(class(segmentsList)[1] != "snpsSegments")
    stop("mismatch segments list, need a list of segments created by the function 'createSegmentsListBySnps' ")
  
  unit <- segmentsList@unit
  step <- segmentsList@gap 
  
  if( unit != "Bases" & unit != "cM")
    stop("Error only cM or Bases are accepted, please make sure to choose between 'cM' or 'Bases' when creating segments list.")
  
    
  
  if(!missing(fileName))
  {
    cat("Warning, you are using an advanced option : a user-specified submap\n")
    res <- readChar(fileName, file.info(fileName)$size)
    res <- unlist(strsplit(res, "\n"))
    submap <- match(res, bedmatrix@snps$id)
  }else{
    submap <- c()
    for(chr in 1:length(segmentsList@snpsSegments))
    {
      map <- segmentsList@snpsSegments[[chr]]
      if(length(map) > 0)
      {
        v <- getMarkerChromosomBySnps(x=bedmatrix, map=map, pas=step, unit=unit) #return an index vector of the marker pick randomly in the segment
        
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

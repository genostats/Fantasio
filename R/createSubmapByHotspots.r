##################################################################################
#This function choose a random marker from every segment of the map element      #
#                                                                                #
#!!! chrSegmentsList : a list of segments of a chromosome                        #
#                                                                                #
#*** return a vector of marker randomly picked                                   #
##################################################################################

getMarkerChromosom <- function(chrSegmentsList)
{
  #find the markers index in the interval
  submap <- c()
  
  #Step 4 : choose a random mkr
  for( i in 1:length(chrSegmentsList))  # throughout the segment
  {
    if(length(chrSegmentsList[[i]]) == 0) next  # empty segment
    if(length(chrSegmentsList[[i]]) == 1) 
    {
      s <- chrSegmentsList[[i]]
    } else { 
      s <- sample(chrSegmentsList[[i]], 1)
    }
    submap <- c(submap, s)
  }
  
  return(submap)
}

#' Creation of a submaps
#' 
#' This function creates a submaps using the list segments created by using the hotspots in the genome
#' 
#' @param bedmatrix a bed.matrix object 
#' @param segmentsList a list of segment for each chromosomes
#' @param epsilon genotype error rate (default is 0.001)
#' @param fileName a fileName with the list of markers wanted for the submap
#' 
#' @details This function will iterates over the list of segments, then for each segments it will pick randomly one marker and put it into a vector.
#' @details Once the iteration over the list is over, the function will create an object and filled some of his slot.
#' @details If you are using a fileName, please make sure to have one marker per line
#' 
#' @return return an hotspots.matrix object with some of his slots filled.
#' 
#' @seealso makeSubmapsByHotspots
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListBySnps(bedMatrix)
#' submaps <- makeSubmapsBySnps(bedMatrix, 5, segmentList)
#' 
#' @export
createSubmapByHotpots <- function(bedmatrix, segmentsList, epsilon = 1e-3, fileName)
{
  if(class(segmentsList)[1] != "hotspot.segments")
    stop("mismatch segments list, need a list of segments created by the function 'createSegmentsListByHotspots' ")
  
  if(!missing(fileName))
  {
    cat("Warning, you are using an advanced option : a user-specified submap\n")
    res <- readChar(fileName, file.info(fileName)$size)
    res <- unlist(strsplit(res, "\n"))
    submap <- match(res, bedmatrix@snps$id)
  }else{
    submap <- c()
    for(chr in 1:length(segmentsList))
    {
      chrMarker <- segmentsList[[chr]]
      randomMarkerVector <- getMarkerChromosom(chrMarker)
      submap <- c(submap, randomMarkerVector)
    }
  }
  
  
  map <- bedmatrix@snps[submap , c("id","chr")]
  if(all(bedmatrix@snps$dist == 0)) {
    map$distance <- bedmatrix@snps$pos[submap]*1e-6
  } else {
    map$distance <- bedmatrix@snps$dist[submap]
  }
  
  log.emiss <- bed.logEmiss(bedmatrix, submap, epsilon)

  new("hotspots.matrix", length(submap), nrow(bedmatrix), submap, 
      bedmatrix@ped[,c("famid", "id", "father", "mother", "sex", "pheno")],
      map, log.emiss, epsilon)

}



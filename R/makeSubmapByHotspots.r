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
  submap <- numeric(length(chrSegmentsList))
  
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
    submap[i] <- s
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
#' @param snpIndices (optional) You have the possibility to pass your own list of markers
#' 
#' @details If snpIndices is given, the function creates a submap corresponding to the given SNPs.
#' @details Otherwise, this function iterates over the list of segments, then for each segments it picks randomly one marker.
#' @details This function is used internally in the package by the function makeAtlasByHotspots
#' 
#' @return return an HostspotsMatrix object with some of his slots filled.
#' 
#' @seealso makeSubmapsByHotspots
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' 
#' @export
makeSubmapByHotspots <- function(bedmatrix, segmentsList, epsilon = 1e-3, snpIndices)
{
  if(class(segmentsList)[1] != "HostspotsSegments")
    stop("mismatch segments list, need a list of segments created by the function 'segmentsListByHotspots' ")
  
  if(!missing(snpIndices))
  {
    submap <- snpIndices
  } else {
    segmentSummary <- segmentsListSummary(segmentsList)

    shift <- cumsum(segmentSummary$number_of_segments)
    shift <- append(0, shift)#0 because I add one after
    max <- shift[length(shift)]

    submap <- numeric(max)

    for(chr in 1:length(segmentsList))
    {
      chrMarker <- segmentsList[[chr]]
      randomMarkerVector <- getMarkerChromosom(chrMarker)
      submap[(shift[chr]+1):shift[chr+1]] <- randomMarkerVector
    }
  }
  
  
  map <- bedmatrix@snps[submap , c("id","chr")]
  if(all(bedmatrix@snps$dist == 0)) {
    map$distance <- bedmatrix@snps$pos[submap]*1e-6
  } else {
    map$distance <- bedmatrix@snps$dist[submap]
  }
  
  log.emiss <- bed.logEmiss(bedmatrix, submap, epsilon)

  new("HostspotsMatrix", length(submap), nrow(bedmatrix), submap, 
      bedmatrix@ped[,c("famid", "id", "father", "mother", "sex", "pheno")],
      map, log.emiss, epsilon)

}



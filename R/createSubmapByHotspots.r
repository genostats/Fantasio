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

##################################################################################
#This function create a submap                                                   #
#                                                                                #
#!!! bedmatrix : a bed.matrix                                                    #                                       
#!!! segmentsList : a list of segments for every chromosomes                     #
#!!! epsilon : the value of epsilon for submaps                                  #
#!!! fileName : a fileName with the list of markers wanted for the submap        #
#                                                                                #
#*** return a new submap object                                                  #
##################################################################################



createSubmapByHotpots <- function(bedmatrix, segmentsList, epsilon = 1e-3, fileName)
{
  
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



#' Creation of a submaps
#' 
#' This function creates a submaps using the list segments created by using gap between markers
#' 
#' @param bedmatrix a bed.matrix object 
#' @param segmentsList a list of segment for each chromosomes
#' @param epsilon the value of epsilon for submaps
#' @param step one marker every step
#' @param unit the unit in which the submaps is to be made, two options are possible "cM" or "Bases"
#' @param fileName You have the possibility to pass your own list of marker, that is to say, your own submaps 
#' 
#' @details This function will iterates over the list of segments, then for each segments it will pick randomly one marker 
#' and from this marker chose a marker every "step" step from back to end.
#' @details Once the iteration over the list is over, the function will create an object and filled some of his slot.
#' @details If you are using a fileName, please make sure to have one marker per line
#' 
#' @return return an snps.matrix object with some of his slots filled.
#' 
#' @seealso \code{\link{makeSubmapsBySnps}}
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListBySnps(bedMatrix)
#' submaps <- makeSubmapsBySnps(bedMatrix, 5, segmentList)
#' 
#' @export
createSubmapBySnps <- function(bedmatrix, segmentsList, epsilon = 1e-3, step=0.5, unit="cM", fileName)
{
  if( unit != "Bases" & unit != "cM")
    stop("Error only cM or Bases are accepted")
  
  if(unit=="Bases")
    step <-  step*1e6
    
  if(segmentsList@gap != step)
      stop("step and gap not equal!")
  
  if(!missing(fileName))
  {
    cat("Warning, you are using an advanced option : a user-specified submap\n")
    res <- readChar(fileName, file.info(fileName)$size)
    res <- unlist(strsplit(res, "\n"))
    submap <- match(res, bedmatrix@snps$id)
  }else{
    submap <- c()
    for(chr in 1:length(segmentsList@snps.segments))
    {
      map <- segmentsList@snps.segments[[chr]]
      if(length(map) > 0)
      {
        v <- getMarkerChromosomBySnps(x=bedmatrix, map=map, pas=step, unit=unit) #return an index vector of the marker pick randomly in the segment
        if(unit == "Bases")
        {
          tmp_dist <- bedmatrix@snps$pos[v]
        }else{
          tmp_dist <- bedmatrix@snps$dist[v]
        }
        #test pour verifier l'ecart entre les snps aux bornes entre deux minis segments
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
  
  new("snps.matrix", step, length(submap), nrow(bedmatrix), submap, 
      bedmatrix@ped[,c("famid", "id", "father", "mother", "sex", "pheno")],
      map, log.emiss, epsilon)
}

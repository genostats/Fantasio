#' Creation of a list of segments 
#' 
#' This function is used to create a list of segments delimited by using hotspots file of the genome.
#' 
#' @param bedmatrix a bed.matrix object 
#' @param intensity Allows the user to select hotspots having recombination intensity higher than a threshold (in cM/Mb). Default value is 10 cM/Mb.
#' @param hotspot_version the version of the hotspots files
#' @param hotspot_file a personnal hotspots file with right columns in it ( see hotspots files )
#' @param verbose Allow the use to know the sequence of action when the list is created, default is TRUE
#' @param number_of_marker a threshold indicating the number of minimum marker in a segment, default is 50
#' 
#' @details This function is used to create a list of chromosome list, which contains segments delimited by using the hotpost file.
#' @details The list is then wrapped under an object of class hotspots.segments for simplicity sake. 
#' @details In each segment you have the index of all the makers in it, these index correpond to markers in the bedmatrix object. 
#' @details The list structure can be analyse using `str()` function on the object (! careful the result can look messy if not handle properly)
#' @details You can give your own hotspots file, be sure to respect the format of the hotpots file, use a default hotpots file for it.
#' 
#' @return an hotspots.segments object
#' 
#' @seealso \code{\link{read.bed.matrix}}
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListByHotspots(bedMatrix)
#' @export
createSegmentsListByHotspots <- function(bedmatrix, intensity = 10 , hotspot_version = "hg19", hotspot_file, verbose = TRUE, number_of_marker = 50)
{
  if(class(bedmatrix)[1] != "bed.matrix")
  {
    stop("Need a bed.matrix to eat")
  }
  
  if(verbose) cat(" Go grab a cup of your favorite beverage, it might takes some time ! :)\n")
  if(verbose) cat(paste("You are currently using version", hotspot_version, "of hotspot\n"))
  
  if(!missing(hotspot_file))
    hotspot <- hotspot_file
  else{
    hotspot <- switch(hotspot_version,
                    hg17 = { data(hotspot_hg17); hotspot_hg17;},
                    hg18 = { data(hotspot_hg18); hotspot_hg18;},
                    hg19 = { data(hotspot_hg19); hotspot_hg19; })
  }
  
  
  #Step 1 : list of all the genome's hotspot
  
  if(verbose) cat("Gathering all hotspots for the genome : ")
  
  segmentsList <- list()
  for ( i in unique(hotspot$Chromosome))
  {
    cat(".")
    chr_hotspot <- hotspot[which(hotspot$Chromosome==i),]
    hotspotsList <- which(chr_hotspot$IntensitycMMb > intensity)
    segment <- cbind(c(0,chr_hotspot$End[hotspotsList]),
                     c(chr_hotspot$Start[hotspotsList],Inf) )
    segmentsList[[i]] <- segment
  }
  if(verbose) cat("\n")
  
  #Step 2 : list of all the marker's position
  
  if(verbose) cat("Gathering all the genome's markers : ")
  
  markersList <- list()
  for( j in unique(hotspot$Chromosome))
  { 
    cat(".")
    markerIndexChromosome <- bedmatrix@snps$pos[bedmatrix@snps$chr==j] 
    markersList[[j]] <- markerIndexChromosome
  }
  cat("\n")
  
  #Step 3 : list of all the segment in the genome

  if(verbose) cat("Finding which markers are between two hotspots : ")
  shift <- sapply(unique(bedmatrix@snps$chr), function(i) which(bedmatrix@snps$chr == i)[1]) - 1L
  
  segmentsAndMarkerList <- list()
  for(i in unique(hotspot$Chromosome))
  {
    cat(".")
    chr_segment <- segmentsList[[i]]
    mkr <- markersList[[i]]
    chr <- list()
    for( j in 1:nrow(chr_segment))
    {
      intercept <- which(mkr > chr_segment[j,1] & mkr < chr_segment[j,2]) #which markers are  between two hotspots
      if (length(intercept)== 0) next
      chr[[j]] <- intercept + shift[[i]]
    }
    segmentsAndMarkerList[[i]] <- chr
    segmentsAndMarkerList[[i]] <- null.remover(segmentsAndMarkerList[[i]])
    segmentsAndMarkerList[[i]] <- cleanHotspots(segmentsAndMarkerList[[i]], number_of_marker)  
  }
  if(verbose) cat("\n")
  
   
  new("hotspot.segments", segmentsAndMarkerList)
} 


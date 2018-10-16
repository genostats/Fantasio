#' Creation of a list of segments 
#' 
#' Creates a list of segments delimited by recombination hotspots
#'
#' @param bedmatrix a bed.matrix object 
#' @param intensity hotspots intensity threshold in cM/Mb
#' @param hotspots a data frame of recombination rates
#' @param minMarkers minimum number of markers in a segment
#' @param verbose if \code{TRUE}, displays information on the process
#' 
#' @details This function creates an object of class hotspots.segments, containing a list of segments delimited
#' by hotspots, The object is a list of list of vectors indices of SNPs. There are as many sublist as
#' chromosomes. The indices correspond to SNPs in \code{bedmatrix}.
#' @details The user can provide a hotspots data frame with a format similar to hotspot_hg19.
#' 
#' @return an hotspots.segments object
#' 
#' @seealso \code{\link{Fantasio}}, \code{\link{segmentsListBySnps}}
#' 
#' @examples  
#' #Please refer to vignette
#'
#' @export
segmentsListByHotspots <- function(bedmatrix, intensity = 10 , hotspots = "hotspot_hg19", minMarkers = 50, verbose = TRUE)
{
  if(class(bedmatrix)[1] != "bed.matrix" )
    stop("Need a bed.matrix")
  
  
  if(verbose) 
    cat( paste("You are currently using version", hotspot_version, "of hotspot\n") )
  
  dataFrameColNames <- c("Chromosome", "Start", "End", "IntensitycMMb")
  if(!all(dataFrameColNames %in% colnames(hotspots)))
    stop("'hotspots' should have columns names  Chromosome, Start, End, IntensitycMMb")
  
  #Step 1 : list of all the genome's hotspot
  
  if(verbose) cat("Gathering all hotspots for the genome : ")
  
  VI <- list()
  for ( i in unique(hotspots$Chromosome))
  {
    cat(".")
    chr_hotspot <- hotspots[which(hotspots$Chromosome==i),]
    w <- which(chr_hotspot$IntensitycMMb > intensity)
    segment <- cbind(c(0,chr_hotspot$End[w]),
                     c(chr_hotspot$Start[w],Inf) )
    VI[[i]] <- segment
  }
  if(verbose) cat("\n")
  
  #Step 2 : list of all the marker's position
  
  if(verbose) cat("Gathering all the genome's markers : ")
  
  VII <- list()
  for(j in unique(hotspots$Chromosome))
  { 
    cat(".")
    v <- bedmatrix@snps$pos[bedmatrix@snps$chr==j] 
    VII[[j]] <- v
  }
  cat("\n")
  
  #Step 3 : list of all the segment in the genome

  if(verbose) cat("Finding which markers are between two hotspots : ")
  shift <- sapply(unique(bedmatrix@snps$chr), function(i) which(bedmatrix@snps$chr == i)[1]) - 1L
  
  VIII <- list()
  for(i in unique(hotspots$Chromosome))
  {
    cat(".")
    chr_segment <- VI[[i]]
    mkr <- VII[[i]]
    chr <- list()
    for( j in 1:nrow(chr_segment))
    {
      b <- which(mkr > chr_segment[j,1] & mkr < chr_segment[j,2]) #which markers are  between two hotspots
      if (length(b)== 0) next
      chr[[j]] <- b + shift[[i]]
    }
    VIII[[i]] <- chr
    VIII[[i]] <- null.remover(VIII[[i]])
    VIII[[i]] <- cleanHotspots(VIII[[i]], number_of_marker)  
  }
  if(verbose) cat("\n")
  
   
  new("HostspotsSegments", VIII)
} 

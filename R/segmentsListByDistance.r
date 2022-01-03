#' Creation of a list of segments 
#' 
#' Creates a list of segments delimited by gaps in the genome.
#' 
#' @param bedmatrix a bed.matrix object 
#' @param gap the minimum gap after which a segment is created
#' @param minMarkers the minimum number of markers in a mini-segment
#' @param nbSegments the number of mini-segments which the segment will be split
#' @param unit the \code{gap} unit, "cM" (centirmorgans) or "Mb" (million of bases)
#' @param verbose if \code{TRUE}, displays information on the process
#' 
#' @details This function is used to create an object of class snpsSegments. It contains a list of segments 
#' delimited by gaps between available markers. Each segment is split in \code{nbSegments} mini-segments.
#' @details A segment is not split if it would result in mini-segments shorter that \code{minMarkers}.
#' The indices in the mini-segments correspond to SNPs in \code{bedmatrix}.
#'
#' @return an snpsSegments object
#' 
#' @seealso \code{\link{Fantasio}}, \code{\link{segmentsListByDistance}}
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' @export
segmentsListByDistance <- function(bedmatrix, gap=0.5, minMarkers=50, nbSegments=20, unit="cM", verbose=TRUE)
{
  if(class(bedmatrix)[1] != "bed.matrix" ) {
    stop("Need a bed.matrix")
  }
  
  if( unit != "Mb" & unit != "cM")
    stop("'unit' should be 'cM' or 'Mb'")

  if(unit == "Mb") {
    gap <- gap * 1e6
  }

  chr.ids <- as.character(intersect( getOption("gaston.autosomes"), unique(bedmatrix@snps$chr)))
  # start and end of a segments
  if(verbose) cat("Finding segments for the genome : ")
  VI <- list()
  for(i in chr.ids) {
    if(verbose) cat(".")
    
    if(unit == "cM") {
      chr_distances <- bedmatrix@snps$dist[which(bedmatrix@snps$chr==i)]
    } else {
      chr_distances <- bedmatrix@snps$pos[which(bedmatrix@snps$chr==i)]
    }

    if(length(chr_distances) == 0) next
    
    k <- c()
    for(j in seq_along(chr_distances))
    {
      if(j == length(chr_distances)) next
      if(chr_distances[j+1] - chr_distances[j] > gap)
        k <- c(k, j, (j+1) )
    }
    
    segment <- cbind( c(0, k[c(FALSE,TRUE)]), 
                      c(k[c(TRUE,FALSE)],Inf))  
    
    
    
    VI[[i]] <- segment
  }
  VI <- VI[!sapply(VI,is.null)]
  
  if(verbose) cat("\n")
  # number of snps in a chr
  VII <- table(bedmatrix@snps$chr)
  VII <- VII[chr.ids]
  
  
  #find the marker of a segment
  if(verbose) cat("Finding which markers are between two segments: ")
  shift <- sapply(chr.ids, function(i) which(bedmatrix@snps$chr == i)[1]) - 1L
  
  VIII <- list()
  for(i in names(VI))
  {
    if(verbose) cat(".")
    chr_segments <- VI[[i]]
    mkr <- seq(1, VII[i])
    chr <- list()
    for(j in seq_len(nrow(chr_segments)))
    {
      b <- which(mkr >= chr_segments[j,1] & mkr <= chr_segments[j,2])
      if(length(b)==0) next
      if(length(b)==1) {
      chr[[j]] <- b + shift[[i]]
      } else {
      c <- c(b[1], b[length(b)])
      chr[[j]] <- c + shift[[i]]
      }
    }
    VIII[[i]] <- chr
  }
  if(verbose) cat("\n")
  
  for(i in seq_along(VIII))
    VIII[[i]] <- null.remover(VIII[[i]])
  
  #finding the mini segments
  if(verbose) cat("Finding mini segments ")
  VIV <- list()
  for(i in seq_along(VIII)) {
    if(verbose) cat(".")
    temp <- list()
    for(j in seq_along(VIII[[i]])) {
    	if (length(VIII[[i]][[j]]) == 1) {
    	temp [[j]] = VIII[[i]][[j]]
    	} 
    	
    	else {
    		if ((length(VIII[[i]][[j]][1]:VIII[[i]][[j]][2]) / nbSegments) >= minMarkers){
        		l <- split(VIII[[i]][[j]][1]:VIII[[i]][[j]][2], ceiling(seq_along(VIII[[i]][[j]][1]:VIII[[i]][[j]][2])/(length(VIII[[i]][[j]][1]:VIII[[i]][[j]][2])/nbSegments))) 
        		for (k in 1:length(l)){
        			if (length(l[[k]]) == 1) next 
        			l[[k]] <- c(l[[k]][1], l[[k]][length(l[[k]])])
        		}
        temp[[j]] <- l
      } else {
        temp[[j]] <- VIII[[i]][[j]]
      }
     }
    }
    VIV[[i]] <- temp
    VIV[[i]] <- null.remover(VIV[[i]])
  }
  if(verbose) cat("\n")
  
  new("snpsSegments", gap, unit, VIV)
}



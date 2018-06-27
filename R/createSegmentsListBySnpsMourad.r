#' Creation of a list of segments 
#' 
#' This function is uses to create a list of segments delimited by using gaps in the genome.
#' 
#' @param bedmatrix a bed.matrix object 
#' @param gap Allows the user to select the minimum gap after which a segment is identified, default is 0.5
#' @param number_of_marker the minimum number of maker in a mini segment, default is 50
#' @param number_of_segments the number of part in which the segment will be splitted
#' @param unit Allow the user to specify in which unit he wants to work
#' @param verbose Allow the use to see the sequance of action going through the algorithm
#' 
#' @details This function is uses to create a list of chromosome list, which contains segments delimited by using the gap between marker.
#' @details The list is then wrap under an object of class snps.segments, which contains two slots, the first slot is the value of the gap used to 
#' create the segment, the second is the list of mini segment and marker in it.
#' @details In each mini segment you have the index of all the makers in it, these index correpond to markers in the bedmatrix object.
#' @details The list structure can be analyze using `str()` function on the object (! careful the result can look messy if not handle properly)
#' @details The number of marker in a mini segment is obtained by computing the length of the segments divided by the number of mini segment wanted. 
#' If the result is less than the number in the variable number_of_marker we keep the segment and do not split it. 
#' Else we use the `split` function to cut the segment into chunks of mini segment. 
#' @details Let us use N as the number of mini segment and T the result of the computation (length of the segment/number of mini segment).
#' If T is an integer the length of the mini segment is exactly T. If T isn't an integer, the length of the mini segment is round(T) or round(T+1).
#' @details The proportion of mini segment with length T+1 increases when T approch round(T+1)
#' @details If (length of the segment/number of mini segment) is equal or greater than the number in the variable number_of_marker, will always have the exact number of 
#' mini segments in the variable number_of_segment.#' 
#' 
#' @return an snps.segments object
#' 
#' @seealso \code{\link{read.bed.matrix}}
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListBySnps(bedMatrix)
#' @export
createSegmentsListBySnps <- function(bedmatrix, gap=0.5, number_of_marker=50, number_of_segments=20, unit="cM", verbose=TRUE,numberOfChromosome=22)
{
  if(class(bedmatrix)[1] != "bed.matrix")
  {
    stop("Need a bed.matrix to eat")
  }
  
  if( unit != "Bases" & unit != "cM")
    stop("Error only cM or Bases are accepted")
  if(unit =="Bases")
    gap <- gap * 1e6
  
  ### start and end of a segments
  if(verbose) cat("Finding segments for the genome : ")
  VI <- list()
  #for(i in unique(bedmatrix@snps$chr)) !!! to be put back
  for(i in 1:numberOfChromosome)
  {
    cat(".")
    chr_distances <- bedmatrix@snps$dist[which(bedmatrix@snps$chr==i)]
    k <- c()
    for(j in 1:length(chr_distances))
    {
      if(j == length(chr_distances))
        next()
      if(chr_distances[j+1] - chr_distances[j] > gap)
        k <- c(k, j, (j+1) )
    }
    segment <- cbind( c(0, k[c(FALSE,TRUE)]), 
                      c(k[c(TRUE,FALSE)],Inf))  
    VI[[i]] <- segment
  }
  
  if(verbose) cat("\n")
  # number of snps in a chr
  VII <- table(bedmatrix@snps$chr)
  
  
  #find the marker of a segment
  if(verbose) cat("Finding which markers are between two segments: ")
  #shift <- sapply(unique(bedmatrix@snps$chr), function(i) which(bedmatrix@snps$chr == i)[1]) - 1L !!!!!
  shift <- sapply(1:numberOfChromosome, function(i) which(bedmatrix@snps$chr == i)[1]) - 1L
  
  VIII <- list()
  #for(i in unique(bedmatrix@snps$chr)) !!!!
  for(i in 1:numberOfChromosome)
  {
    cat(".")
    chr_segments <- VI[[i]]
    mkr <- seq(1, VII[i])
    chr <- list()
    for(j in 1:nrow(chr_segments))
    {
      b <- which(mkr >= chr_segments[j,1] & mkr <= chr_segments[j,2])
      if(length(b)==0) next()
      chr[[j]] <- b + shift[[i]]
    }
    VIII[[i]] <- chr
  }
  if(verbose) cat("\n")
  
  for(i in 1:length(VIII))
    VIII[[i]] <- null.remover(VIII[[i]])
  
  #finding the mini segments
  if(verbose) cat("Finding mini segments ")
  VIV <- list()
  for(i in 1:length(VIII))
  {
    cat(".")
    temp <- list()
    for(j in 1:length(VIII[[i]]))
    {
      if((length(VIII[[i]][[j]]) / number_of_segments) >= number_of_marker ) #>= number_of_marker in one segments 
      {
        #decouper le segment en N (=number_of_segments) mini-segments de taille T (=length(VIII[[i]][[j]])/number_of_segments) marqueurs.
        #si T est un entier, la taille du mini-segment est exactement T
        #si T n'est pas un entier, la taille est round(T) ou round(T)+1 en fonction des resultats de ceiling(...)
        #La proportion de mini-segments de taille round(T)+1 augmente lorsque T approche de round(T)+1
        #Note: cette commande donne toujours exactement N (=number_of_segments) mini-segments
        l <- split(VIII[[i]][[j]], ceiling(seq_along(VIII[[i]][[j]])/(length(VIII[[i]][[j]])/number_of_segments))) 
        temp[[j]] <- l
      }else{
        temp[[j]] <- VIII[[i]][[j]]
      }
    }
    VIV[[i]] <- temp
    VIV[[i]] <- null.remover(VIV[[i]])
  }
  if(verbose) cat("\n")
  
  new("snps.segments", gap, VIV)
}





#Enlever car n'est plus necessaire suite a la ligne de commande precedente. A confirmer?
#if(length(l) > 1)
#{
#  index <- c()
#  for(k in 1:length(l))
#  {
#    if(length(l[[k]]) <= number_of_marker)
#    {
#      if(k == length(l))
#        l[[k-1]] <- c(l[[k]], l[[k-1]])
#      else
#        l[[k+1]] <- c(l[[k]], l[[k+1]])
#        #l[[k]] <- NULL
#    }else{
#      index <- c(index, k)
#    }
#  }
#  l <- l[index]
#}

createSegmentsListBySnpsIsuru <- function(bedmatrix, gap=0.5, number_of_marker=50, unit="cM", verbose=TRUE)
{
  if( unit != "Bases" & unit != "cM")
    stop("Error only cM or Bases are accepted")
  if(unit =="Bases")
    gap <- pas * 1000
  
  ### start et end d'un segment
  if(verbose) cat("Finding segments for the genome : ")
  VI <- list()
  for(i in 1:22)
  {
    cat(".")
    chr <- bedmatrix@snps$dist[which(bedmatrix@snps$chr==i)]
    c <- c()
    for(j in 1:length(chr))
    {
      if(j == length(chr))
        next()
      if(chr[j+1] - chr[j] > gap)
        c <- c(c, j, (j+1) )
    }
    segment <- cbind(c(0, c[c(F,T)]), c(c[c(T,F)],Inf))  
    VI[[i]] <- segment
  }
  if(verbose) cat("\n")
  #distance pour chaque marqueur
  if(verbose) cat("Gathering all the genome's markers : ")
  VII <- list()
  for(i in 1:22)
  {
    cat(".")
    VII[[i]] <- bedmatrix@snps$dist[bedmatrix@snps$chr==i]
  }
  if(verbose) cat("\n")
  #trouver les marqueurs dans un segment
  if(verbose) cat("Finding which markers are between two segments: ")
  shift <- sapply(1:22, function(i) which(bedmatrix@snps == i)[1]) - 1L
  VIII <- list()
  for(i in 1:22)
  {
    cat(".")
    chr_segments <- VI[[i]]
    mkr <- VII[[i]]
    chr <- list()
    for(j in 1:nrow(chr_segments))
    {
      b <- which(mkr > chr_segments[j,1] & mkr < chr_segments[j,2])
      if(length(b)==0)next()
      chr[[j]] <- b + shift[[i]]
    }
    VIII[[i]] <- chr
  }
  if(verbose) cat("\n")
  #supprimer les NULL
  for(i in 1:length(VIII))
    VIII[[i]] <- null.remover(VIII[[i]])

  #trouver les minis-segments
  if(verbose) cat("Finding mini segments ")
  VIV <- list()
  for(i in 1:length(VIII))#parcourt de tous les chromosomes
  {
    cat(".")
    temp <- list()
    for(j in 1:length(VIII[[i]]))#parcourt de tous les segmets des chromosomes
    {
      #split(VIII[[i]][[j]], ceiling(seq_along(d)/20))
      l <- split(VIII[[i]][[j]], ceiling(seq_along(VIII[[i]][[j]])/number_of_marker))##besoin d'un vecteur en premier arg
      #l <- split(VIII[[i]][[j]], cut(seq_along(VIII[[i]][[j]]), number_of_segments, labels = FALSE)) 
      if(length(l) > 1)
      {
        for(k in 1:length(l))
        {
          if(length(l[[k]]) < number_of_marker)
          {
            l[[k-1]] <- c(l[[k-1]], l[[k]])
            l[[k]] <- NULL
            
          }
        }
      }
      temp[[j]] <- l
    }
    VIV[[i]] <- temp
  }
  if(verbose) cat("\n")
  
  new("hotspot.segments", VIV)
}
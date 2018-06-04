createSegmentsListBySnps <- function(bedmatrix, gap=0.5, number_of_marker=50, number_of_segments=20, unit="cM", verbose=TRUE)
{
  if( unit != "Bases" & unit != "cM")
    stop("Error only cM or Bases are accepted")
  if(unit =="Bases")
    gap <- gap * 1e6
  
  ### start et end d'un segment
  if(verbose) cat("Finding segments for the genome : ")
  VI <- list()
  #for(i in unique(bedmatrix@snps$chr))
  for(i in 1:22)
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
  # nbre de SNPs par chr
  VII <- table(bedmatrix@snps$chr)
  
  
  #trouver les marqueurs dans un segment
  if(verbose) cat("Finding which markers are between two segments: ")
  #shift <- sapply(unique(bedmatrix@snps$chr), function(i) which(bedmatrix@snps$chr == i)[1]) - 1L
  shift <- sapply(1:22, function(i) which(bedmatrix@snps$chr == i)[1]) - 1L
  
  VIII <- list()
  #for(i in unique(bedmatrix@snps$chr))
  for(i in 1:22)
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
      #length of the segment / number_of_segments wanted >= number_of_marker in one segments 
      if((length(VIII[[i]][[j]]) / number_of_segments) >= number_of_marker )
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
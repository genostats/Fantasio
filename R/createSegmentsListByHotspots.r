# x = an bedmatrix, n = the number of submap

createSegmentsListByHotspots <- function(bedmatrix, intensity = 10 , hotspot_version = "hg19", hotspot_file, verbose = TRUE, number_of_marker = 50)
{
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
  ##on boucle seulement sur le nombre de chromosome dans le fichier de hotspots, si plus de chromosome dans les donnees => infos perdues
  if(verbose) cat("Gathering all hotspots for the genome : ")
  
  VI <- list()
  for ( i in unique(hotspot$Chromosome))
  {
    cat(".")
    chr_hotspot <- hotspot[which(hotspot$Chromosome==i),]
    w <- which(chr_hotspot$IntensitycMMb > intensity)
    segment <- cbind(c(0,chr_hotspot$End[w]),
                     c(chr_hotspot$Start[w],Inf) )
    VI[[i]] <- segment
  }
  if(verbose) cat("\n")
  
  #Step 2 : list of all the marker's position
  
  if(verbose) cat("Gathering all the genome's markers : ")
  
  VII <- list()
  for( j in unique(hotspot$Chromosome))
  { 
    cat(".")
    v <- bedmatrix@snps$pos[bedmatrix@snps$chr==j] 
    VII[[j]] <- v
  }
  cat("\n")
  
  #Step 3 : list of all the segment in the genome
  #TODO renommer les variables !!!
  if(verbose) cat("Finding which markers are between two hotspots : ")
  shift <- sapply(unique(bedmatrix@snps$chr), function(i) which(bedmatrix@snps$chr == i)[1]) - 1L
  
  VIII <- list()
  for(i in unique(hotspot$Chromosome))
  {
    cat(".")
    chr_segment <- VI[[i]]
    mkr <- VII[[i]]
    chr <- list()
    for( j in 1:nrow(chr_segment))
    {
      b <- which(mkr > chr_segment[j,1] & mkr < chr_segment[j,2])#which markers are  between two hotspots
      if (length(b)== 0) next
      chr[[j]] <- b + shift[[i]]
    }
    VIII[[i]] <- chr
    VIII[[i]] <- null.remover(VIII[[i]])
    VIII[[i]] <- cleanHotspots(VIII[[i]], number_of_marker)  
  }
  if(verbose) cat("\n")
  
   
  new("hotspot.segments", VIII)
} 


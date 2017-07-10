# x = an bedmatrix, n = the number of submap

segments <- function(x, intensity = 10 , hotspot_version = "hg19", verbose = TRUE)
{
  if(verbose) cat(" Go grab a cup of your favorite beverage, it might takes some time ! :)\n")
  if(verbose) cat(paste("You are currently using version", hotspot_version, "of hotspot\n"))
  hotspot <- switch(hotspot_version,
                    hg17 = { data(hotspot_hg17); hotspot_hg17;},
                    hg18 = { data(hotspot_hg18); hotspot_hg18;},
                    hg19 = { data(hotspot_hg19); hotspot_hg19; })
  
  
  #Step 1 : list of all the genome's hotspot
  
  if(verbose) cat("Gathering all hotspots for the genome : ")
  
  VI <- list()
  for ( i in 1:22)
  {
    cat(".")
    chr_hotspot <- hotspot[which(hotspot[,1]==i),]
    w <- which(chr_hotspot$IntensitycMMb > intensity)
    segment <- cbind(c(0,chr_hotspot$End[w]),
                     c(chr_hotspot$Start[w],Inf) )
    VI[[i]] <- segment
  }
  if(verbose) cat("\n")
  
  #Step 2 : list of all the genome's markers
  
  if(verbose) cat("Gathering all the genome's markers : ")
  
  VII <- list()
  for( j in 1:22)
  { 
    cat(".")
    v <- x@snps$pos[x@snps$chr==j] #lapply(x@snps$pos[x@snps$chr==j], submap.error)
    VII[[j]] <- v
  }
  cat("\n")
  
  #Step 3 : list of all the segment in the genome
  #TODO renommer les variables !!!
  if(verbose) cat("Finding which markers are between two hotspots : ")
  shift <- sapply(1:22, function(i) which(x@snps == i)[1]) - 1L
  
  VIII <- list()
  for(i in 1:22)
  {
    cat(".")
    chr_segment <- VI[[i]]
    mkr <- VII[[i]]
    chr <- list()
    for( j in 1:nrow(chr_segment))
    {
      b <- which(mkr > chr_segment[j,1] & mkr < chr_segment[j,2])#which markers are contains between two hotspots
      if (length(b)== 0) next
      chr[[j]] <- b + shift[[i]]
    }
    VIII[[i]] <- chr
  }
  if(verbose) cat("\n")
 
  new("hotspot.segments",VIII)
} 


# x = an bedmatrix, n = the number of submap

submap <- function(x, n = 100, intensity = 10 , hotspot_version = "hg19")
{
  cat(" Go grab a cup of your favorite beverage, it might takes some time ! :)")
  cat("\n")
  cat(paste("You are currently using version", hotspot_version, "of hotspot"))
  cat("\n")
  hotspot <- switch(hotspot_version,
                    hg17 = { data(hotspot_hg17); hotspot_hg17;},
                    hg18 = { data(hotspot_hg18); hotspot_hg18;},
                    hg19 = { data(hotspot_hg19); hotspot_hg19; })
  
  
  #Step 1 : list of all the genome's hotspot
  
  cat("Gathering all informations about genome's hotspot : ")
  
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
  cat("\n")
  
  #Step 2 : list of all the genome's markers
  
  cat("Gathering all the genome's markers : ")
  
  VII <- list()
  for( j in 1:22)
  { 
    cat(".")
    v <- x@snps$pos[x@snps$chr==j] #lapply(x@snps$pos[x@snps$chr==j], submap.error)
    VII[[j]] <- v
  }
  cat("\n")
  
  #Step 3 : list of all the segment in the genome
  
  cat("Gathering all the segments for the genome thanks to previous infos : ")
  shift <- sapply(1:22, function(i) which(x@snps == i)[1]) - 1
  
  VIII <- list()
  for(i in 1:22)
  {
    cat(".")
    chr_segment <- VI[[i]]
    mkr <- VII[[i]]
    chr <- list()
    for( j in 1:nrow(chr_segment))
    {
      b <- which(mkr > chr_segment[j,1] & mkr < chr_segment[j,2])
      if (length(b)== 0) next
      
      chr[[j]] <- b + shift[[i]]
      
    }
    VIII[[i]] <- chr
  }
  cat("\n")
  
  submap <- array(list(), c(n,1))
  for ( i in 1:n)
   { 
     set.seed(i)  # FOR DEBUG PURPOSE TO BE REMOVED!!!
     cat("creating submap number : ", i,"/", n, "\n" )
     spider <- createSubmap(x, VIII) #invisible(createSubmap(x, VIII))
     submap[[i,1]] <- spider
   }
 #  # #for ( i in 1:n)
 # #{ 
 #   set.seed(92)  # FOR DEBUG PURPOSE TO BE REMOVED!!!
 #   cat("creating submap number : ", 92,"/", n, "\n" )
 #   spider <- createSubmap(x, VIII) #invisible(createSubmap(x, VIII))
 #   submap[[i,1]] <- spider
 # #}
  cat("Creation of all the Submap over ! \n")
  return(submap)
  
  
  #x@Submap <- Submap 
  #return(x)
}


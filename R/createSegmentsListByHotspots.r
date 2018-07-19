#' Creation of a list of segments 
#' 
#' This function is used to create a list of segments delimited by using hotspots files of the genome.
#' 
#' @param bedmatrix a bed.matrix object 
#' @param intensity Allows the user to select hotspots having recombination intensity higher than a threshold in cM/Mb (default is 10)
#' @param hotspot_version the version of the hotspots files (default is "hg19")
#' @param hotspot_file a personnal hotspots file with right columns in it ( see hotspots files )
#' @param verbose whether you want informations about the computation process (default is TRUE)
#' @param number_of_marker a threshold indicating the number of minimum marker in a segment (default is 50)
#' 
#' @details This function is used to create a list of chromosome list, which contains segments delimited by using the hotpost file.
#' @details The list is then wrapped under an object of class hotspots.segments for simplicity sake. 
#' @details In each segment you have the index of all the makers in it, these index correpond to markers in the bedmatrix object. 
#' @details The list structure can be analyse using `str()` function on the object (! careful the result can look messy if not handle properly)
#' @details You can give your own hotspots file, be sure to respect the format of the hotpots file, use a default hotpots file for it.
#' 
#' @return an hotspots.segments object
#' 
#' @seealso read.bed.matrix
#' @seealso Fantasio
#' 
#' @examples  
#' #install.packages("HGDP.CEPH", repos="https://genostats.github.io/R/") ## make this only one time
#' require(Fantasio)
#' require(HGDP.CEPH)
#' filepath <-system.file("extdata", "hgdp_ceph.bed", package="HGDP.CEPH")
#' x <- read.bed.matrix(filepath)
#' x <- set.stats(x)
#' x.me <- select.inds(x, population == "Bedouin")
#' x.me@ped$pheno <- rep(2,48)
#' segmentList <- createSegmentsListByHotspots(x.me)
#' @export
createSegmentsListByHotspots <- function(bedmatrix, intensity = 10 , hotspot_version = "hg19", hotspot_file, verbose = TRUE, number_of_marker = 0)
{
  if(class(bedmatrix)[1] != "bed.matrix" )
    stop("Need a bed.matrix to eat")
  
  
  if(verbose) 
    cat( paste("You are currently using version", hotspot_version, "of hotspot\n") )
  
  if(!missing(hotspot_file))
  {
    dataFrameColNames <- c("Chromosome", "Start", "End", "IntensitycMMb")
    test <- sapply(dataFrameColNames, function(i) i %in% colnames(hotspot_file))
    
    if(all(test))
      hotspot <- hotspot_file
    else
      stop("Your file is not written properly, please ensure that you have a dataframe with atleast the following columns names : 
            Chromosome, Start, End, IntensitycMMb and the correct informations in it 
            Thank you.")
  }else{
    hotspot <- switch(hotspot_version,
                    hg17 = Fantasio::hotspot_hg17,
                    hg18 = Fantasio::hotspot_hg18,
                    hg19 = Fantasio::hotspot_hg19 )
  }
  
  
  #Step 1 : list of all the genome's hotspot
  
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
      b <- which(mkr > chr_segment[j,1] & mkr < chr_segment[j,2]) #which markers are  between two hotspots
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


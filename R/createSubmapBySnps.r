##################################################################################
#This function create a submaps using a step for creating a submap                #
#                                                                                #
#                                                                                #
#*** return a new submap object                                                  #
##################################################################################
#map : segment et mini_segment du chr 

createSubmapBySnps <- function(bedmatrix, segmentsList, epsilon = 1e-3, step=0.5, unit="cM", fileName)
{
  if( unit != "Bases" & unit != "cM")
    stop("Error only cM or Bases are accepted")
  
  if(unit=="Bases")
    step <-  step*1e6
    
  if(segmentsList@gap != step)
      stop("step and gap not equal!")
  
  if(!missing(fileName))
  {
    cat("Warning, you are using an advanced option : a user-specified submap\n")
    res <- readChar(fileName, file.info(fileName)$size)
    res <- unlist(strsplit(res, "\n"))
    submap <- match(res, bedmatrix@snps$id)
  }else{
    submap <- c()
    for(chr in 1:length(segmentsList@snps.segments))
    {
      map <- segmentsList@snps.segments[[chr]]
      if(length(map) > 0)
      {
        v <- getMarkerChromosomBySnps(bedmatrix, map, step, unit)
        
        if(unit == "Bases")
        {
          tmp_dist <- bedmatrix@snps$pos[v]
        }
        else{
          tmp_dist <- bedmatrix@snps$dist[v]
        }
        #test pour verifier l'ecart entre les snps aux bornes entre deux mini segments
        tmp <- diff(tmp_dist)
        v <- v[which(tmp >= step)]
        submap <- c(submap, v)
      }
    }
  }
  
  
  map <- bedmatrix@snps[submap , c("id","chr")]
  if(unit == "cM")
  {
    if(all(bedmatrix@snps$dist == 0)) {
      map$distance <- bedmatrix@snps$pos[submap]*1e-6
    } else {
      map$distance <- bedmatrix@snps$dist[submap]
    }
  }
  else
  {
    if(all(bedmatrix@snps$pos == 0)) {
      map$distance <- bedmatrix@snps$dist[submap]*1e6
    } else {
      map$distance <- bedmatrix@snps$pos[submap]
    }
  }
  
  
  log.emiss <- bed.logEmiss(bedmatrix, map=submap, epsilon=epsilon)
  
  new("snps.matrix", step, length(submap), nrow(bedmatrix), submap, 
      bedmatrix@ped[,c("famid", "id", "father", "mother", "sex", "pheno")],
      map, log.emiss, epsilon)
}

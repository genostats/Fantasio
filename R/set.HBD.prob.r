#x an matrix with a and f computed 
set.HBD.prob <- function(submaps, list.id, quality = 95)
{
  if(!missing(list.id))
  {
    if(list.id == "all")
    {
      condition <- 1:nrow(submaps@submap_summary)
    }else{
      vec <- strsplit(list.id, "_")
      condition <- sapply(vec, function(i) which(submaps@submap_summary$FID == i[1] & submaps@submap_summary$IID == i[2]))
    }
  }else{
    condition <- which(submaps@submap_summary$QUALITY >= quality & submaps@submap_summary$INBRED)
    if(length(condition) == 0)
    {
      cat("WARNING :No inbred found with the following default parameters : QUALITY = ", quality , "; inbred = TRUE 
          you can try to change quality parameters or give a vector of individual or use 'all' parameter \n 
          !!! Instead using all the individuals in the sample\n")
      #for simplicity sake, not returning an empty results
      condition <- 1:nrow(submaps@submap_summary)
    }
      
  }
  
  id    <- as.vector(submaps@submap_summary$IID[condition])
  famid <- as.vector(submaps@submap_summary$FID[condition])
  
  for(i in 1:length(submaps@atlas))
  {
    HBD_prob <- matrix(NA, nrow = length(condition), ncol = submaps@atlas[[i]]@ncol)#HBD matrix
    dimnames(HBD_prob) <- list(rownames(HBD_prob) <- paste(famid,id, sep = "_"), colnames(HBD_prob) <- c(submaps@atlas[[i]]@map$id))
    
    for (j in 1:nrow(HBD_prob))
    {
      j1 <- condition[j]
      if(!is.na(submaps@atlas[[i]]@a[j1]) & (submaps@atlas[[i]]@a[j1] <= 1) & !is.na(submaps@atlas[[i]]@f[j1]))
      {
        HBD_prob[j,1:ncol(HBD_prob)] <-forward.backward(get.log.emiss(submaps@atlas[[i]], j1), 
                                                        submaps@atlas[[i]]@delta.dist, 
                                                        submaps@atlas[[i]]@a[j1],
                                                        submaps@atlas[[i]]@f[j1] )[2,]
      }
    }
    submaps@atlas[[i]]@HBD.prob <- HBD_prob 
  }
  l <- list(submaps, condition)
  return(l)
}





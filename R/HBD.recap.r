
HBD.recap <- function(h)
{
  snps <- as.vector(sapply(h, function(x) colnames(x@HBD.prob))) #every markers that has been choosen on submaps
  b <- as.data.frame(table(snps),  stringsAsFactors=FALSE)  #count the number of times the marker has been choosen
  n <- matrix(nrow = nrow(b), ncol = (h[[1]]@nrow) + 1) # create the matrix : HBD.recap
  namevector <- c("Freq", rownames(h[[1]]@HBD.prob)) #vector or column names
  dimnames(n) <- list(rownames(n) <- c(b$snps), colnames(n) <- namevector) #give name to matrix
  n[,1] <- b$Freq
  a <- ncol(n) #use to put data on the matrix
  
  
  az <- vapply(h, function(hh) match( b$snps, colnames(hh@HBD.prob)), integer(length(b$snps))  )
  cat("This is gonna takes some time pls get yourself comfortable :) \n")
  #through all the markers
  for( i in 1:nrow(n))
  {
    az1 <- az[i,]
    I <- which(!is.na(az1))
    #v <- vapply( I, function(i) h[[i]]@HBD.prob[, az1[i]], double(nrow(h[[1]]@HBD.prob)) )
    v <- sapply( I, function(i) h[[i]]@HBD.prob[, az1[i]]) 
    n[i,seq(2,a, by = 1)] <- apply(v,1,mean)
  }  
  return(n)
}  

#for( i in 1:nrow(n))
#{
#  cat(i,  ",")
#  id <- b$snps[i]
#  v <- Reduce( cbind, lapply( h, function(hh) if(id %in% colnames(hh@HBD.prob)) hh@HBD.prob[,str()] ))
#  n[i,seq(2,a, by = 1)] <- apply(v,1,mean)
#}

#v <- vapply( I, function(i) h[[i]]@HBD.prob[, az1[i]], double(nrow(h[[i]]@HBD.prob)) )
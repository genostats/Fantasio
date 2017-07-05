
essai <- function(h)
{
  snps <- as.vector(sapply(h, function(x) colnames(x@HBD.prob)))
  b <- as.data.frame(table(snps),  stringsAsFactors=FALSE)
  
  n <- matrix(nrow = nrow(b), ncol = (h[[1]]@nrow) + 1)
  namevector <- c("Freq", rownames(h[[1]]@HBD.prob))
  dimnames(n) <- list(rownames(n) <- c(b$snps), colnames(n) <- namevector)
  
  n[,1] <- b$Freq
  
  #through all the markers
  for( i in 1:nrow(n))
  {
    id <- b$snps[i]
    v <- Reduce( cbind, lapply( h, function(hh) if(id %in% colnames(hh@HBD.prob)) hh@HBD.prob[,id] ))
    
    
    res <- c()
    
    #through all individuals
    for( j in 1:nrow(v))
    {
      res <- c(res, mean(v[j,]))
    }
    
    n[i,seq(2,a, by = 1)] <- res
  }
  
  
  
  
}  

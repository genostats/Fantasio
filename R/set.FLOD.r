set.FLOD <- function(submaps, condition, q = 1e-4)
{
  #Computation of FLOD score with the formula
  for(i in 1:length(submaps@atlas))
  {
    FLOD_prob <- matrix(0.0, nrow = nrow(submaps@atlas[[i]]@HBD.prob), ncol = submaps@atlas[[i]]@ncol)
    dimnames(FLOD_prob) <- list(rownames(FLOD_prob) <- rownames(submaps@atlas[[i]]@HBD.prob), colnames(FLOD_prob) <- colnames(submaps@atlas[[i]]@HBD.prob))
    
    for (j in 1:nrow(FLOD_prob))
    {
      if( (submaps@atlas[[i]]@a[condition[j]] < 1) & !is.na(submaps@atlas[[i]]@f[condition[j]]) )
      {
        FLOD_prob[j,] <- log10((submaps@atlas[[i]]@HBD.prob[j,] + q * ( 1 - submaps@atlas[[i]]@HBD.prob[j,]))/
                                 (submaps@atlas[[i]]@f[condition[j]] + q * ( 1 - submaps@atlas[[i]]@f[condition[j]]))) 
      }
    }
    submaps@atlas[[i]]@FLOD <- FLOD_prob 
  }
  return(submaps)
  
}

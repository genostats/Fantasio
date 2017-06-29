submap.summary <- function(h, a.threshold = 1)
{  
  f <-  sapply(h, function(x) x@f) 
  a <-  sapply(h, function(x) x@a)
  p <- sapply(h , function(x) x@p.lrt)
  f[a > a.threshold] <- NA
  a[a > a.threshold] <- NA
  p[a > a.threshold] <- NA
  pLRT_MEDIAN <- apply(p, 1, median, na.rm=TRUE)
  
  l <- p < 0.05
  
  nValidSubmap <- c()
  for (i in 1:nrow(l))
  {
    nValidSubmap <- c(nValidSubmap, sum(l[i,], na.rm = TRUE ))
  }
  
  
  df <- data.frame(FID           = h[[1]]@ped$famid, 
                   IID           = h[[1]]@ped$id,
                   STATUS        = h[[1]]@ped$pheno,
                   SUBMAPS       = length(h),
                   QUALITY       = rowSums(sapply(h, function(x) x@a < a.threshold) )/length(h),
                   F_MIN         = apply(f, 1, min, na.rm = TRUE), 
                   F_MAX         = apply(f, 1, max, na.rm = TRUE),
                   F_MEAN        = apply(f, 1, mean, na.rm = TRUE),
                   F_MEDIAN      = apply(f, 1, median, na.rm=TRUE),
                   A_MEDIAN      = apply(a, 1, median, na.rm=TRUE),
                   pLRT_MEDIAN   = pLRT_MEDIAN, 
                   INBRED        = pLRT_MEDIAN < 0.05, 
                   LRT_inf_0.05  = nValidSubmap)
  
  #write.table(df,"dataframe.txt",sep="\t",col.names = colnames(df))
  #library("xlsx", character.only = TRUE)
  write.xlsx(df, "dataframe.xlsx") 
  return(df)
  
  
}


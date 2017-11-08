submap.summary <- function(submaps, a.threshold = 1)
{  
  f <-  sapply(submaps, function(x) x@f) 
  a <-  sapply(submaps, function(x) x@a)
  p <- sapply(submaps , function(x) x@p.lrt)
  w.a <- (a > a.threshold)
  f[w.a] <- NA
  p[w.a] <- NA
  a[w.a] <- NA
  pLRT_MEDIAN <- apply(p, 1, median, na.rm=TRUE)
  
  l <- p < 0.05
  
  nValidSubmap <- c()
  for (i in 1:nrow(l))
  {
    nValidSubmap <- c(nValidSubmap, sum(l[i,], na.rm = TRUE ))
  }
  
  #treat the case when quality is equal to NA
  quality <- (rowSums(sapply(submaps, function(x) x@a < a.threshold) )*100)/length(submaps)
  n <- which(is.na(quality))
  quality[n] <- 0
  
  
  df <- data.frame(FID           = submaps[[1]]@ped$famid, 
                   IID           = submaps[[1]]@ped$id,
                   STATUS        = submaps[[1]]@ped$pheno,
                   SUBMAPS       = paste(length(submaps), "/", length(submaps)),
                   QUALITY       = quality,
                   F_MIN         = apply(f, 1, min, na.rm = TRUE), 
                   F_MAX         = apply(f, 1, max, na.rm = TRUE),
                   F_MEAN        = apply(f, 1, mean, na.rm = TRUE),
                   F_MEDIAN      = apply(f, 1, median, na.rm=TRUE),
                   A_MEDIAN      = apply(a, 1, median, na.rm=TRUE),
                   pLRT_MEDIAN   = pLRT_MEDIAN, 
                   INBRED        = pLRT_MEDIAN < 0.05, 
                   pLRT_inf_0.05 = nValidSubmap)
  
  for(i in 1:nrow(df))
  {
    if(!is.finite(df$F_MEDIAN[i]))
      df[i,5:13] <- NA
  }
  
  #write.table(df,"dataframe.txt",sep="\t",col.names = colnames(df))
  #library("xlsx", character.only = TRUE)
  #write.xlsx(df, "dataframe.xlsx") 
  return(df)
  
  
}


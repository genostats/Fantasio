#' Submap summary
#' 
#' This function creates a summary on the submaps created. 
#' 
#' @param submaps a list of submaps 
#' @param a.threshold  the maximum value for a (default is 1)
#' 
#' @details This function gives for each genotyped individual summary statistics about the calculations.
#' @details This function returns a dataframe with 13 columns :
#' @details - FID : family identifier
#' @details - IID : individual identifier
#' @details - STATUS : status (1 non affected, 2 affected, 0 unknown)
#' @details - submaps : number of submaps used
#' @details - QUALITY: percentage of valid submaps (i.e. submaps with a < 1)
#' @details - F_MIN: minimum f on valid submaps
#' @details - F_MAX: maximum f on valid submaps
#' @details - F_MEAN: mean f on valid submaps
#' @details - F_MEDIAN: median f on valid submaps (recommended to estimate f)
#' @details - A_MEDIAN: median a on valid submaps (recommended to estimate a)
#' @details - pLRT_MEDIAN: median p-value of LRT tests on valid submaps
#' @details - INBRED: a flag indicating if the individual is inbred (pLRT_MEDIAN <0.05) or not
#' @details - pLRT_<0.05: number of valid submaps with a LRT having a p-value below 0.05

#' @return this function returns a dataframe.
#' 
#' @seealso setSummary
#' 
#' @export
submapSummary <- function(submaps, a.threshold = 1)
{  
  if(class(submaps[[1]])[1] != "snps.matrix" & class(submaps[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.") 
  
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
  submaps_used <- rowSums(sapply(submaps, function(x) x@a <= a.threshold))
  quality <- (submaps_used*100)/length(submaps)
  n <- which(is.na(quality))
  quality[n] <- 0
  
  
  df <- data.frame(FID           = submaps[[1]]@ped$famid, 
                   IID           = submaps[[1]]@ped$id,
                   STATUS        = submaps[[1]]@ped$pheno,
                   SUBMAPS       = paste(submaps_used, "/", length(submaps)),
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
  
  return(df)
}


segment.summary <- function(s)
{
  #number of segments
  n_seg <- c()
  for(i in 1:length(s))
  {
    n_seg <- c(n_seg, length(s[[i]]))
  }
  n_seg[23] <- NA
  n_seg <- n_seg[!is.na(n_seg)]
  
  
  #number of markers 
  n_mark <- c()
  for(i in 1:length(s))
  {
    res <- c()
    for(j in 1:length(s[[i]]))
    {
      res <- c(res, length(s[[i]][[j]]))
    }
    res <- sum(res)
    n_mark <- c(n_mark, res)
  }
  
  
  #dataframe
  df <- data.frame(
    chromosom = 1:22,
    number_of_segments = n_seg, 
    number_of_markers= n_mark
  )
  df
}


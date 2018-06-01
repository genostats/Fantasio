summary.marker <- function(submaps, bedmatrix)
{
  b <- summary.map(submaps)
  res <- c()
  for(i in 1:length(submaps))
  {
    res <- c(res, length(which(b$Freq == i)))
  }
  zero <- length(bedmatrix@snps$chr) - sum(res)
  
  df <- data.frame(
    number_of_time_picked = 0:length(submaps),
    number_of_markers = c(zero, res)
  )
  df
}
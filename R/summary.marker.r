summary.marker <- function(h)
{
  b <- summary.map(h)
  res <- c()
  for(i in 1:length(h))
  {
    res <- c(res, length(which(b$Freq == i)))
  }
  zero <- length(x@snps$chr) - sum(res)
  
  df <- data.frame(
    number_of_time_picked = 0:length(h),
    number_of_markers = c(zero, res)
  )
  df
}
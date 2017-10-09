summary.marker <- function(h, bedmatrix)
{
  b <- summary.map(h@atlas)
  res <- c()
  for(i in 1:length(h@atlas))
  {
    res <- c(res, length(which(b$Freq == i)))
  }
  zero <- length(x@snps$chr) - sum(res)
  
  df <- data.frame(
    number_of_time_picked = 0:length(h@atlas),
    number_of_markers = c(zero, res)
  )
  df
}
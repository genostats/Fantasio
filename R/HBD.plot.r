#plot for one chromosom and one individual
HBD.plot.chr <- function(x, chr, ind) 
{
  #find the distance
  w <- which(x@map$chr == chr)
  d <- x@map$distance[w]
  n <- length(d)
  r <- x@HBD.prob[ind,w]
  
  
  #plotting
  ymax <- max(1, max(x@HBD.prob[ind, w], na.rm = TRUE)) 
  xmax <- max(d, na.rm = TRUE)
  ymin <- 0
  #d=distance for each chr i ; r=HFLOD for each individual j & chr i
  plot(d, r, type="l", pch=16, xlim = c(0,xmax), ylim = c(ymin, ymax), 
       xlab = "Position (cM)", ylab="HBD", main= paste("HBD (individual ",ind,")", sep=""))  
  
  abline(h=seq(0, 1, 0.1), col="grey", lwd=1.5, lty=2)
}
# paralleliser !!
submaps <- function(x, n = 100, segments, epsilon = 1e-3, run.festim = TRUE, verbose = TRUE) {

  submap <- vector("list", n)
  for ( i in 1:n)
  { 
    set.seed(i)  # FOR DEBUG PURPOSE TO BE REMOVED!!!
    if(verbose) cat("creating submap number : ", i,"/", n, "\n" )
    spider <- createSubmap(x, segments, epsilon) 
    if(run.festim) 
      spider <- festim(spider)
    submap[[i]] <- spider
  }

  if(verbose) cat("Creation of all the Submaps over ! \n")
  return(submap)
}



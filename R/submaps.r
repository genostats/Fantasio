
submaps <- function(x, n = 100, segment = segments(x), n.cores = 1, epsilon = 1e-3, run.festim = TRUE, proba = TRUE,  verbose = TRUE) {

  ff <- function(i, run.festim) {
    set.seed(i)  # FOR DEBUG PURPOSE TO BE REMOVED!!!
    spider <- createSubmap(x, segment, epsilon) 
    if(run.festim) 
      spider <- festim(spider, probs = proba, verbose = verbose)
    spider
  }

  if(n.cores != 1 & .Platform$OS.type != "unix") {
    warning("FORK cluster unavailable only one core used")
    n.cores <- 1
  }
  
  if(n.cores == 1) {
    submap <- lapply(1:n, ff, run.festim = run.festim)
  }else {
    RNGkind("L'Ecuyer-CMRG")
    s <- matrix(.Random.seed, nrow = 1)
    for(i in 2:n.cores) 
      s <- rbind(s, nextRNGStream(s[i-1,]))
    cl <- makeForkCluster(n.cores) #create slaves
    parLapply(cl, 1:n.cores, function(i) .Random.seed <<- s[i,] ) #use of paralelisation
    submap <- parLapply(cl, 1:n, ff, run.festim = run.festim)
    stopCluster(cl)
    gc()
  }
  
  #submap <- vector("list", n)
  #for ( i in 1:n)
  #{ 
  #  set.seed(i)  # FOR DEBUG PURPOSE TO BE REMOVED!!!
  #  if(verbose) cat("creating submap number : ", i,"/", n, "\n" )
  #  spider <- createSubmap(x, segment, epsilon) 
  #  if(run.festim) 
  #    spider <- festim(spider)
  #  submap[[i]] <- spider
  #}

  if(verbose) cat("Creation of all the Submaps over ! \n")
  submaps <- new("list.submaps", submap, x)
  submaps <- setSummary(submaps)
  submaps
}



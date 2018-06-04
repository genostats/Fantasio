makeSubmapsBySnps <- function(bedmatrix, n = 100, segmentsList = createSegmentsListBySnps(x), n.cores = 1, epsilon = 1e-3, run.festim=TRUE, list.id, proba=TRUE,
                    recap.by.segments = FALSE,  verbose=TRUE, debug=FALSE, step=0.5, unit="cM", threshold=0.5, q = 1e-4, quality=95, n.consecutive.marker=5) {
  
  if(class(segmentsList)[1] != "snps.segments")
    stop("mismatch segments list, need a list of segments created by the function 'createSegmentsListBySnps' ")
  
  ff <- function(i, run.festim) {
    spider <- createSubmapBySnps(bedmatrix, segmentsList, epsilon, step, unit) 
    if(run.festim) 
      spider <- festim(spider, verbose=verbose, debug=debug)
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

  if(missing(list.id))
  {
    
  }
  
  submaps <- new("list.submaps", submap, bedmatrix, segmentsList@snps.segments, recap.by.segments)
  submaps <- setSummary(submaps, run_a_f = run.festim, probs = proba, by_segments=recap.by.segments, list.id=list.id, threshold=threshold, q=q, quality=quality, n.consecutive.marker=n.consecutive.marker)
  if(verbose) cat("Creation of all the Submaps over ! \n")
  submaps
}

##################################################################################
#This function creates n submaps and allows the creation of summary files        #
#                                                                                #
#!!! bedmatrix : a bed.matrix                                                    #
#!!! n : the number of submaps wanted                                            #             
#!!! segmentsList : a list of segment for each chromosomes                       #
#!!! n.cores : the number of cores to use if you want to compute submaps using   #
#    parellelism                                                                 #
#!!! epsilon : the value of epsilon for submaps                                  #
#!!! run.festim: whether you want to computes a, f, p.lrt, likelihood0/1 for     #
#    each submaps                                                                #        
#!!! list.id : you can either                                                    #
#             - ignore this parameter if you want to compute HBD, FLOD and HFLOD #
#                  for individuals who are considerated INBRED and with a QUALITY#
#                  greater or equal to 95%                                       #
#              - enter a list of individual for a computation of HBD, FLOD score #
#                 HFLOD score for them                                           #
#              - the character "all" for a computation of HBD, FLOD score and    #
#                  HFLOD score for every individual                              #
#!!! run.proba : whether you want to computes HBD, FLOD score and HFLOD score    #
#!!! recap.by.segments : if you want the summary of probabilities by snps or by  #
#                        by segments                                             #
#!!! verbose : whether you want informations about computations                  #
#                                                                                #
#*** return a new list object containing every dataframe and object created      #
##################################################################################

makeSubmapsByHotspots <- function(bedmatrix, n = 100, segmentsList = createSegmentsListByHotspots(bedmatrix), n.cores = 1, epsilon = 1e-3, run.festim=TRUE, list.id, run.proba=TRUE,
                    recap.by.segments = FALSE,  verbose=TRUE, debug=FALSE, threshold=0.5, q = 1e-4, quality=95, n.consecutive.marker=5) {

  if(class(segmentsList)[1] != "hotspot.segments")
    stop("mismatch segments list, need a list of segments created by the function 'createSegmentsListByHotspots' ")
  
  ff <- function(i, run.festim) {
    spider <- createSubmapByHotpots(bedmatrix, segmentsList, epsilon=epsilon) 
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
  
  submaps <- new("list.submaps", submap, bedmatrix, segmentsList, recap.by.segments)
  submaps <- setSummary(submaps, run_a_f = run.festim, probs = run.proba, by_segments=recap.by.segments, list.id=list.id, threshold=threshold, q=q, quality=quality, n.consecutive.marker=n.consecutive.marker)
  if(verbose) cat("Creation of all the Submaps over ! \n")
  submaps
}



#' Creation of submaps using gap between markers
#' 
#' This function creates n submaps and allows the creation of summary files
#' 
#' @param bedmatrix a bed.matrix object 
#' @param n the number of submaps wanted(default is 100)
#' @param segmentsList a list of segment for each chromosomes
#' @param n.cores the number of cores to use if you want to compute submaps using parellelism (default is 1)
#' @param epsilon genotype error rate (default is 0.001)
#' @param run.festim whether you want to computes a, f, p.lrt, likelihood0/1 for each submaps (default is TRUE)
#' @param list.id a list of individuals (see details for more informations)
#' @param run.proba whether you want to computes HBD, FLOD score and HFLOD score (default is TRUE)  
#' @param recap.by.segments if you want the summary of probabilities by snps or by segments (default is FALSE)
#' @param verbose whether you want informations about computations (default is TRUE)
#' @param debug whether you want advanced output for the computation process (default is FALSE)
#' @param threshold the value of the threshold when finding HBD segment, threshold is the probability of being HBD or not (default is 0.5)
#' @param q Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param quality Allows the user to choose the minimal quality (in \%) to include an inbred individual into the analysis (default is 95)
#' @param n.consecutive.marker the number of consecutive marker with a probabilitie equal or greater to the value of threshold, to be use to fing HBD segments
#' 
#' 
#' @details This function is used to create submaps by randomly picking 
#' a marker in each segment in the list of segments given by the argument segmentsList.
#' @details After the creation of one submap is finished, the function pass 
#' it to `festim`, a function that will computes different values (a, f, p.lrt, likelihood0, likelihood1)
#' @details After that we then creates the different summary to interpret 
#' easily the results of the computation. 
#' @details This is done by the function `setSummary`. This function will 
#' call all of the function neccessary to obtained the different summary needed for our data.
#' @details You can check all the summary after they are been created by accesing their slots.
#' @details When using recap.by.segments with true value, we then consider 
#' that the snps picked randomly in a segment is a representant of that segment.
#' @details you can pass several arguments to the list.id arguments : 
#'\itemize{
#'  \item{ignore this parameter if you want to compute HBD, FLOD and HFLOD for INRED individuals and with a QUALITY}
#'  \item{enter a list of individual for a computation of HBD, FLOD 
#'  score HFLOD score for them, each element of the vector should contains the familyId_individualId}
#'  \item{the character "all" for a computation of HBD, FLOD score and HFLOD score for every individual}
#' }
#' @details When doing a number N of submaps using the options recap.by.segments, 
#' @details the different values computed from the 
#' @details snps picked randomly (a, f, p.lrt, ...) in a segments is considered different values 
#' @details for the same segments.
#' 
#' @return return a new list object containing every dataframe and object created 
#' 
#' @seealso Fantasio
#' @seealso makeSubmapsBySnps
#' @seealso createSegmentsListByHotspots
#' @seealso festim
#' @seealso set.HBD.prob
#' @seealso set.FLOD
#' @seealso set.HFLOD
#' @seealso recap
#' @seealso setSummary
#' @seealso submapLikelihood
#' @seealso submapEstim
#' @seealso summaryMarker
#' @seealso submapSummary
#' @seealso HBD.segments
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' 
#' @export
makeAllSubmapsBySnps <- function(bedmatrix, n = 100, segmentsList = createSegmentsListBySnps(bedmatrix), n.cores = 1, epsilon = 1e-3,
                                 run.festim=TRUE, list.id, run.proba=TRUE, recap.by.segments = FALSE, verbose=TRUE, debug=FALSE,
                                 threshold=0.5, q = 1e-04, quality=95, n.consecutive.marker=5)
{
  if(class(segmentsList)[1] != "snps.segments")
    stop("mismatch segments list, need a list of segments created by the function 'createSegmentsListBySnps' ")
  
  ff <- function(i, run.festim) {
    spider <- createSubmapBySnps(bedmatrix, segmentsList, epsilon=epsilon) 
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
    cl <- makeForkCluster(n.cores) #creation of slaves
    parLapply(cl, 1:n.cores, function(i) .Random.seed <<- s[i,] ) #use of paralelisation
    submap <- parLapply(cl, 1:n, ff, run.festim = run.festim)
    stopCluster(cl)
    gc()
  }

  submaps <- new("list.submaps", submap, bedmatrix, segmentsList@snps.segments, recap.by.segments, segmentsList@unit , segmentsList@gap)
  
  submaps <- setSummary(submaps, run_a_f = run.festim, probs = run.proba, recap.by.segments=recap.by.segments, list.id=list.id, threshold=threshold, q=q, quality=quality, n.consecutive.marker=n.consecutive.marker)
  if(verbose) cat("Creation of all the Submaps by distance over ! \n")
  submaps
}

#' Creation of submaps based on hotposts in the genome
#' 
#' This function creates N submaps and allows the creation of summary files
#' 
#' @param bedmatrix a bed.matrix object 
#' @param n the number of submaps wanted(default is 100)
#' @param segmentsList a list of segment for each chromosomes
#' @param n.cores the number of cores to use if you want to compute submaps using parellelism (default is 1)
#' @param epsilon genotype error rate (default is 0.001)
#' 
#' 
#' @details This function is used to create submaps by randomly picking a marker
#'  in each segment in the list of segments given by the argument segmentsList.
#' @details After the creation of one submap is finished, the function pass it 
#' to `festim`, a function that will computes different values (a, f, p.lrt, likelihood0, likelihood1)
#' @details After that we then creates the different summary to interpret easily 
#' the results of the computation. 
#' @details This is done by the function `setSummary`. This function will call all 
#' of the function neccessary to obtained the different summary needed for our data.
#' @details You can check all the summary after they are been created by accesing their slots.
#' @details When using recap.by.segments with true value, we then consider that the 
#' snps picked randomly in a segment is a representant of that segment.
#' @details When doing a number N of submaps using the options recap.by.segments, 
#' @details the different values computed from the 
#' @details snps picked randomly (a, f, p.lrt, ...) in a segments is considered different values 
#' @details for the same segments.
#' 
#' @return return a new list object containing every dataframe and object created 
#' 
#' @seealso Fantasio
#' @seealso makeAtlasByDistance
#' @seealso segmentsListByHotspots
#' @seealso festim
#' @seealso setHBDprob
#' @seealso setFLOD
#' @seealso setHFLOD
#' @seealso recap
#' @seealso setSummary
#' @seealso submapLikelihood
#' @seealso submapEstim
#' @seealso submapSummary
#' @seealso HBDsegments
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' @export
makeAtlasByHotspots <- function(bedmatrix, n = 100, segmentsList = segmentsListByHotspots(bedmatrix), n.cores = 1, epsilon = 1e-3) {

  if(class(segmentsList)[1] != "HostspotsSegments")
    stop("mismatch segments list, need a list of segments created by the function 'segmentsListByHotspots' ")
  
  ff <- function(i) {
    makeSubmapByHotspots(bedmatrix, segmentsList, epsilon = epsilon) 
  }

  if(n.cores != 1 & .Platform$OS.type != "unix") {
    warning("FORK cluster unavailable only one core used")
    n.cores <- 1
  }
  
  if(n.cores == 1) {
    submap <- lapply(seq_len(n), ff)
  } else {
    RNGkind("L'Ecuyer-CMRG")
    s <- matrix(.Random.seed, nrow = 1)
    for(i in 2:n.cores) 
      s <- rbind(s, nextRNGStream(s[i-1,]))
    cl <- makeForkCluster(n.cores) 
    parLapply(cl, seq_len(n.cores), function(i) .Random.seed <<- s[i,] ) 
    submap <- parLapply(cl, seq_len(n), ff)
    stopCluster(cl)
    gc()
  }
  
  new("atlas", submap, bedmatrix, segmentsList, NA)
}

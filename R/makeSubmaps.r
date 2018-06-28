#' Creation of submaps using hotposts in the genome
#' 
#' This function creates N submaps and allows the creation of summary files
#' 
#' @param bedmatrix a bed.matrix object 
#' @param n the number of submaps wanted  
#' @param segmentsList a list of segment for each chromosomes
#' @param n.cores the number of cores to use if you want to compute submaps using parellelism 
#' @param epsilon the value of epsilon for submaps
#' @param run.festim whether you want to computes a, f, p.lrt, likelihood0/1 for each submaps  
#' @param list.id you can either :
#'     - ignore this parameter if you want to compute HBD, FLOD and HFLOD 
#'       for individuals who are considerated INBRED and with a QUALITY
#'       greater or equal to 95%}
#'     - enter a list of individual for a computation of HBD, FLOD score HFLOD score for them
#'     - he character "all" for a computation of HBD, FLOD score and HFLOD score for every individual
#' @param run.proba whether you want to computes HBD, FLOD score and HFLOD score  
#' @param recap.by.segments if you want the summary of probabilities by snps or by segments
#' @param verbose whether you want informations about computations
#' @param debug whether you want advanced output for the computation process
#' @param threshold the value of the threshold when finding HBD segment
#' @param q Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual. Default is 0.0001.
#' @param quality Allows the user to choose the minimal quality (in %) to include an inbred individual into the analysis. Default is 95 (%).
#' @param n.consecutive.marker the number of consecutive marker with a probabilitie equal or greater to the value of threshold, to be use to fing HBD segments
#' @param byHotspots whether the submaps are to be made using segments created by hotspots or segments created by gap between markers
#' @param step the value of the step to be made to chose the next marker in each segment
#' @param unit this argument is use only when you make submaps using segments created using the gap between markers, two options are allowed "Base", "cM"
#' 
#' @details byHotspots == TRUE
#' @details This function is used to create submaps by randomly picking a marker in each segment in the list of segments given by the argument segmentsList.
#' @details After the creation of one submap is finished, the function pass it to `festim`, a function that will computes different values (a, f, p.lrt, likelihood0, likelihood1)
#' @details After that we then creates the different summary to interpret easily the results of the computation. 
#' @details This is done by the function `setSummary`. This function will call all of the function neccessary to obtained the different summary needed for our data.
#' @details You can check all the summary after they are been created by accesing their slots.
#' @details When using recap.by.segments with true value, we then consider that the snps picked randomly in a segment is a representant of that segment.
#' @details When doing a number N of submaps the different values computed from the snps picked randomly (a, f, p.lrt, ...) in a marker is considered different values 
#' @details for the same segments.
#' 
#' @return return a new list object containing every dataframe and object created 
#' 
#' @seealso \code{\link{createSegmentsListByHotspots}}
#' @seealso \code{\link{createSegmentsListBySnps}}
#' @seealso \code{\link{festim}}
#' @seealso \code{\link{set.HBD.prob}}
#' @seealso \code{\link{set.FLOD}}
#' @seealso \code{\link{set.HFLOD}}
#' @seealso \code{\link{recap}}
#' @seealso \code{\link{setSummary}}
#' @seealso \code{\link{submapLikelihood}}
#' @seealso \code{\link{submapEstim}}
#' @seealso \code{\link{summaryMarker}}
#' @seealso \code{\link{submapSummary}}
#' @seealso \code{\link{HBD.segments}}
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListByHotspots(bedMatrix)
#' submaps <- makeSubmap(bedMatrix, 5, segmentList)
#' @export
makeSubmaps <- function(bedmatrix, n = 100, segmentsList, n.cores = 1, epsilon = 1e-3, run.festim=TRUE, list.id, run.proba=TRUE,
                    recap.by.segments = FALSE,  verbose=TRUE, debug=FALSE, threshold=0.5, q = 1e-4, quality=95, n.consecutive.marker=5, byHotspots=TRUE,
                    step = 0.5, unit = "cM") {

  statusTwo <- which(bedmatrix@ped$pheno == 2)
  if(length(statusTwo) == 0)
  {
    stop("Error no individual with a STATUS equal to 2. ")
  }
  
  if(byHotspots)
  {
    if(class(segmentsList)[1] != "hotspot.segments")
    stop("mismatch segments list, need a list of segments created by the function 'createSegmentsListByHotspots' ")
  }else{
    if(class(segmentsList)[1] != "snps.segments")
      stop("mismatch segments list, need a list of segments created by the function 'createSegmentsListBySnps' ")
  }
  
  ff <- function(i, run.festim, byHotspots) {
    if(byHotspots)
      oneSubmap <- createSubmapByHotpots(bedmatrix, segmentsList, epsilon=epsilon) 
    else 
      oneSubmap <- createSubmapBySnps(bedmatrix, segmentsList, epsilon=epsilon, step=step, unit=unit) 
    if(run.festim) 
      oneSubmap <- festim(oneSubmap, verbose=verbose, debug=debug)
    oneSubmap
  }

   if(n.cores != 1 & .Platform$OS.type != "unix") {
     warning("FORK cluster unavailable only one core used")
     n.cores <- 1
   }
   
   if(n.cores == 1) {
     submap <- lapply(1:n, ff, run.festim = run.festim, byHotspots = byHotspots)
   }else {
     RNGkind("L'Ecuyer-CMRG")
     s <- matrix(.Random.seed, nrow = 1)
     for(i in 2:n.cores) 
       s <- rbind(s, nextRNGStream(s[i-1,]))
     cl <- makeForkCluster(n.cores) #create slaves
     parLapply(cl, 1:n.cores, function(i) .Random.seed <<- s[i,] ) #use of paralelisation
     submap <- parLapply(cl, 1:n, ff, run.festim = run.festim, byHotspots = byHotspots)
     stopCluster(cl)
     gc()
   }
  
  

  if(byHotspots)
    submaps <- new("list.submaps", submap, bedmatrix, segmentsList, recap.by.segments)
  else 
    submaps <- new("list.submaps", submap, bedmatrix, segmentsList@snps.segments, recap.by.segments)
  submaps <- setSummary(submaps, run_a_f = run.festim, probs = run.proba, by_segments=recap.by.segments, list.id=list.id, threshold=threshold, q=q, quality=quality, n.consecutive.marker=n.consecutive.marker, test=statusTwo)
  if(verbose) cat("Creation of all the Submaps over ! \n")
  submaps
}



#' Wrapper for the package Fantasio
#' 
#' This fonctions is used as a wrapper for the package Fantasio to create segments list and submaps
#' 
#' @param bedmatrix a bed.matrix object 
#' @param segments The method to use for submap creation ("Hotspots" or "Distance")
#' @param segment.options a  list of arguments to the function that will create the segments list. The names attribute of args gives the argument names. 
#' @param n the number of submaps wanted  
#' @param n.cores the number of cores to use if you want to compute submaps using parellelism 
#' @param epsilon the value of epsilon for submaps
#' @param run.festim whether you want to computes a, f, p.lrt, likelihood0/1 for each submaps  
#' @param list.id a list of individuals from which to computes HBD, 
#' FLOD score and HFLOD score (see details for more informations) 
#' @param run.proba whether you want to computes HBD, FLOD score and HFLOD score  
#' @param recap.by.segments if you want the summary of probabilities by snps or by segments
#' @param verbose whether you want informations about computations
#' @param debug whether you want advanced output for the computation process
#' @param threshold the value of the threshold when finding HBD segment
#' @param q Allows the user to choose the assumed frequency of the mutation
#'  involved in the disease for each individual. Default is 0.0001.
#' @param quality Allows the user to choose the minimal quality 
#' to include an inbred individual into the analysis. Default is 95.
#' @param n.consecutive.marker the number of consecutive marker with a 
#' probabilities equal or greater to the value of threshold, to be use to fing HBD segments
#' 
#' 
#' @details This function is a wrapper to make the usage of Fantasio, the package, easier. The function calls two differents function : 
#' @details `createSegmentsListBySnps` and `createSegmentsListBySnps`. The first function `createSegmentsListBySnps` is used to create 
#' @details a list of segments though the genome. The second function `createSegmentsListBySnps` is used to create submaps.
#' @details The segments arguments accept only two options : Hotspots or Distance.
#' @details Depending of the value of the segments argument (Hotspots or Distance) the segments are made by using the hotspots of recombinaison
#' @details in the genome thanks to hotspots file, then the submaps are made by Hotspots by picking randomly a maker in every segments created before.
#' @details Or the segments are made by using the gaps between markers, then the submaps are made by picking a random marker in
#' @details every segments and going through each segment from left to right using a supply step (by default it is 0.5 cM).
#' @details To see which arguments can be passed in the list of args in the segment.options please go check the createSegmentListBySnps/Hotspots functions manual.
#' 
#' @seealso read.bed.matrix
#' @seealso createSegmentsListBySnps
#' @seealso makeAllSubmapsBySnps
#' @seealso makeAllSubmapsByHotspots
#' @seealso makeAllSubmapsByHotspots
#' 
#' @export
Fantasio <- function(bedmatrix, segments = c("Hotspots", "Distance"), segment.options, 
                     n = 100, n.cores = 1, epsilon = 0.001, run.festim = TRUE, list.id, 
                     run.proba = TRUE, recap.by.segments = FALSE, verbose = TRUE, 
                     debug = FALSE, threshold = 0.5, q = 1e-04, quality = 95, 
                     n.consecutive.marker = 5) {
  
  
  qdfuohfuis3489dsifs93iksdi2_tmp <<- bedmatrix #to correct : nodes produced errors; first error: external pointer is not valid
  segments <- match.arg(segments)
  if(missing(segment.options)) # default value
    segment.options <- list()
  
  if(segments == "Distance" & recap.by.segments) {
    recap.by.segments <- FALSE
    warning('segments = "Distance" implies recap.by.segments = FALSE')
  }
  
  if(segments == "Distance") {
    s <- do.call(createSegmentsListBySnps, c(bedmatrix = bedmatrix, segment.options))
    h <- makeAllSubmapsBySnps(qdfuohfuis3489dsifs93iksdi2_tmp, n, s, n.cores, epsilon, run.festim, list.id, run.proba, recap.by.segments, verbose, debug, threshold, q, quality, n.consecutive.marker)
  } else {
    s <- do.call(createSegmentsListByHotspots, c(bedmatrix = bedmatrix, segment.options))
    h <- makeAllSubmapsByHotspots(qdfuohfuis3489dsifs93iksdi2_tmp, n, s, n.cores, epsilon, run.festim, list.id, run.proba, recap.by.segments, verbose, debug, threshold, q, quality, n.consecutive.marker)
  }
  
  rm(qdfuohfuis3489dsifs93iksdi2_tmp, pos = ".GlobalEnv")
  h
}
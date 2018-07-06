#' Creation of submaps using gap between markers
#' 
#' Creates n submaps and generates summaries following the specified options
#' 
#' @param bedmatrix a bed.matrix object 
#' @param segments The method to use for submap creation ("Hotspots" or "Distance")
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
#' @details manuel à finir !!!!!

Fantasio <- function(bedmatrix, segments = c("Hotspots", "Distance"), segment.options, 
    n = 100, n.cores = 1, epsilon = 0.001, run.festim = TRUE, list.id, 
    run.proba = TRUE, recap.by.segments = FALSE, verbose = TRUE, 
    debug = FALSE, threshold = 0.5, q = 1e-04, quality = 95, 
    n.consecutive.marker = 5) {

  segments <- match.arg(segments)
  if(missing(segment.options)) # valeurs par défauts
    segment.options <- list()

  if(segments == "Distance" & recap.by.segments) {
    recap.by.segments <- FALSE
    warning('segments = "Distance" implies recap.by.segments = FALSE')
  }

  if(segments == "Distance") {
    s <- do.call(createSegmentsListBySnps, c(bedmatrix = bedmatrix, segment.options))
    h <- makeAllSubmapsBySnps(x, n, s, n.cores, epsilon, run.festim, list.id, run.proba, recap.by.segments, verbose, debug, threshold, q, quality, n.consecutive.marker)
  } else {
    s <- do.call(createSegmentsListByHotspots, c(bedmatrix = bedmatrix, segment.options))
    h <- makeAllSubmapsByHotspots(x, n, s, n.cores, epsilon, run.festim, list.id, run.proba, recap.by.segments, verbose, debug, threshold, q, quality, n.consecutive.marker)
  }
  h
}

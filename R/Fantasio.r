#' Wrapper for the package Fantasio
#' 
#' This function is used as a wrapper for the package Fantasio to create segments list and submaps
#' 
#' @param bedmatrix a bed.matrix object 
#' @param segments The method to use for submap creation ("Hotspots" or "Distance")
#' @param segment.options a  list of arguments to the function that will create the segments list. The names attribute of args gives the argument names. 
#' @param n the number of submaps wanted(default is 100)
#' @param n.cores the number of cores to use if you want to compute submaps using parellelism (default is 1)
#' @param epsilon genotype error rate (default is 0.001)
#' @param run.proba whether you want to computes HBD, FLOD score and HFLOD score (default is TRUE)  
#' @param recap.by.segments if you want the summary of probabilities by snps or by segments (default is FALSE)
#' @param verbose whether you want informations about computations (default is TRUE)
#' @param threshold the value of the threshold when finding HBD segment, threshold is the probability of being HBD or not (default is 0.5)
#' @param q Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param quality Allows the user to choose the minimal quality (in \%) to include an inbred individual into the analysis (default is 95)
#' @param n.consecutive.markers the number of consecutive marker with a probabilitie equal or greater to the value of threshold, to be use to fing HBDsegments
#' 
#' 
#' @details This function is a wrapper to make the usage of Fantasio, the package, easier. The function calls differents function : 
#' @details The first function `segmentsListByDistance` and`segmentsListByHotspots` is used to create a list of segments though the genome. 
#' @details The second function `makeAtlasBySnsp`and `makeAtlasByHotspots` is used to create submaps.
#' @details The segments arguments accept only two options : Hotspots or Distance.
#' @details Depending of the value of the segments argument (Hotspots or Distance) the segments are made by using the hotspots of recombinaison
#' @details in the genome thanks to hotspots file, then the submaps are made by Hotspots by picking randomly a maker in every segments created before.
#' @details Or the segments are made by using the gaps between markers (option Distance), then the submaps are made by picking a random marker in
#' @details every segments and going through each segment from left to right using a supply step (by default it is 0.5 cM).
#' @details To see which arguments can be passed in the list of args in the segment.options please go check the `segmentsListByDistance` or `segmentsListByHotspots` functions manual.
#' @details When doing a number N of submaps using the options recap.by.segments, 
#' @details the different values computed from the 
#' @details snps picked randomly (a, f, p.lrt, ...) in a segments is considered different values 
#' @details for the same segments.
#' 
#' @seealso read.bed.matrix
#' @seealso segmentsListByDistance
#' @seealso makeAtlasByDistance
#' @seealso makeAtlasByHotspots
#' @seealso makeAtlasByHotspots
#'
#' @examples
#' #Please refer to vignette 
#'
#' 
#' @export
Fantasio <- function (bedmatrix, segments = c("Hotspots", "Distance"), segment.options,
                      n = 100, n.cores = 1, epsilon = 0.001,
                      run.proba = TRUE, recap.by.segments = FALSE, verbose = TRUE,
                      threshold = 0.5, q = 1e-04, quality = 95,
                      n.consecutive.markers = 5)
{
  segments <- match.arg(segments)
  if (missing(segment.options))
    segment.options <- list()
  if(!("verbose" %in% segment.options))
    segment.options$verbose = FALSE

  if (segments == "Distance" & recap.by.segments) {
    recap.by.segments <- FALSE
    warning("segments = \"Distance\" implies recap.by.segments = FALSE")
  }
  if (segments == "Distance") {
    s <- do.call(segmentsListByDistance, c(bedmatrix = bedmatrix, segment.options))
    h <- makeAtlasByDistance(get(deparse(substitute(bedmatrix))), n, s, n.cores, epsilon)
    h <- festim(h, n.cores = n.cores, verbose = verbose)
    h <- setSummary(h, probs = run.proba, recap.by.segments = recap.by.segments, threshold = threshold, 
                    q = q, quality = quality, n.consecutive.markers = n.consecutive.markers)
  }
  else {
    s <- do.call(segmentsListByHotspots, c(bedmatrix = bedmatrix, segment.options))
    h <- makeAtlasByHotspots(get(deparse(substitute(bedmatrix))), n, s, n.cores, epsilon)
    h <- festim(h, n.cores = n.cores, verbose = verbose)
    h <- setSummary(h, probs = run.proba, recap.by.segments = recap.by.segments, threshold = threshold, 
                    q = q, quality = quality, n.consecutive.markers = n.consecutive.markers)
  }
  h
}

#' Creation summary statistic files
#' 
#' This function is uses to ouput all the summary files after the creation of the atlas by the following two functions : 
#' `makeAtlasByHotspots` and `makeAtlasByDistance` is over.
#' 
#' 
#' @param atlas a atlas object
#' @param list.id you can either :
#'     - ignore this parameter if you want to compute HBD, FLOD and HFLOD 
#'       for individuals who are considerated INBRED and with a QUALITY
#'       greater or equal to 95%}
#'     - enter a list of individual for a computation of HBD, FLOD score HFLOD score for them
#'     - use "all" for a computation of HBD, FLOD score and HFLOD score for every individual
#' @param probs a flag indicating if the HBD probabilities and FLOD score has been computed for the submaps (default is TRUE)
#' @param recap.by.segments if the summary files has to be computed considering segments or snps (defaut is FALSE) 
#' (more information in the documentation of `makeAtlasByHotspots` or`makeAtlasByDistance` functions)
#' @param threshold the value of the threshold when finding HBD segment, threshold is the probability of being HBD or not (default is 0.5)
#' @param q  Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param quality The minimum percentage use to assume if a submap is valid (default is 95)
#' @param n.consecutive.markers the number of consecutive markers with a probability equal or greater to the value of the threshold, to be used to find HBDsegments (default is 5)
#' 
#' @details the function filled the empty slots of the atlas object. It calls different functions and uses the results of each one to filled 
#' the object. This function computes summary for the following elements :
#' @details - summary for likelihood0/likelihood1
#' @details - summary for estimation of a and f
#' @details - marker summary
#' @details - summary for the submaps
#' @details - summary for HBD probabilities
#' @details - summary for HFLOD
#' @details - comutation of HBDsegments
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
#' 
#' @export
setSummary <- function (atlas, list.id, probs = TRUE, recap.by.segments = FALSE,
    q = 1e-04, threshold = 0.5, quality = 95, n.consecutive.markers = 5) {

  if(class(atlas)[1] != "atlas")
    stop("Need an atlas")
 
  atlas@likelihood_summary <- submapLikelihood(atlas@submaps_list)
  atlas@estimation_summary <- submapEstim(atlas@submaps_list)
  atlas@marker_summary <- summaryMarker(atlas)
  atlas@submap_summary <- suppressWarnings(submapSummary(atlas@submaps_list))
  
  test <- which(atlas@bedmatrix@ped$pheno == 2)
  if(probs & length(test) == 0 & missing(list.id) ) {
    warning("No cases (pheno = 2). HBD, FLOD and HFLOD are computed on all individuals")
    list.id <- "all"
  }
  
  if(probs) {
    l1 <- setHBDprob(atlas, list.id=list.id, quality=quality)
    atlas <- l1[[1]]
    atlas <- setFLOD(submaps=atlas, condition=l1[[2]], q=q)
   
    # test class of first submap to check if recap is ok 
    if(class(atlas@submaps_list[[1]])[1] == "snpsMatrix" & recap.by.segments) {
      warning("Submaps created by Distance force 'recap.by.segments = FALSE'")
      recap.by.segments <- FALSE
    }
    atlas@bySegments <- recap.by.segments
    l2 <- recap(atlas, list.id = list.id)
    atlas@HBD_recap <- l2[[1]]
    atlas@FLOD_recap <- l2[[2]]
    atlas@HBDsegments <- HBDsegments(atlas, threshold = threshold, n.consecutive.markers = n.consecutive.markers) 
    atlas@HFLOD <- setHFLOD(atlas)
  }
  atlas
}


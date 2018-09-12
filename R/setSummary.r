#' Creation summary statistic files
#' 
#' This function is uses to ouput all the summary files after the creation of the submaps by the following two functions : 
#' `makeSubmapsByHotspots` and `makeSubmapsBySnps` is over.
#' 
#' 
#' @param submaps a submapsList object
#' @param list.id you can either :
#'     - ignore this parameter if you want to compute HBD, FLOD and HFLOD 
#'       for individuals who are considerated INBRED and with a QUALITY
#'       greater or equal to 95%}
#'     - enter a list of individual for a computation of HBD, FLOD score HFLOD score for them
#'     - use "all" for a computation of HBD, FLOD score and HFLOD score for every individual
#' @param run_a_f a flag indicating if the estimation of a and f has been computed for the submaps using `festim` (default is TRUE)
#' @param probs a flag indicating if the HBD probabilities and FLOD score has been computed for the submaps (default is TRUE)
#' @param recap.by.segments if the summary files has to be computed considering segments or snps (defaut is FALSE) 
#' (more information in the documentation of `makeSubmapsByHotspots` or`makeSubmapsBySnps` functions)
#' @param threshold the value of the threshold when finding HBD segment, threshold is the probability of being HBD or not (default is 0.5)
#' @param q  Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param quality The minimum percentage use to assume if a submap is valid (default is 95)
#' @param n.consecutive.marker the number of consecutive marker with a probabilitie equal or greater to the value of the threshold, to be use to fing HBDsegments (default is 5)
#' 
#' @details the function filled the empty slots of the submapsList object. It calls different functions and uses the results of each one to filled 
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
#' @seealso makeSubmapsBySnps
#' @seealso createSegmentsListByHotspots
#' @seealso festim
#' @seealso setHBDprob
#' @seealso setFLOD
#' @seealso setHFLOD
#' @seealso recap
#' @seealso setSummary
#' @seealso submapLikelihood
#' @seealso submapEstim
#' @seealso summaryMarker
#' @seealso submapSummary
#' @seealso HBDsegments
#' 
#' 
#' @export
setSummary <- function (submaps, list.id, run_a_f = TRUE, probs = TRUE, recap.by.segments = FALSE,
    q = 1e-04, threshold = 0.5, quality = 95, n.consecutive.marker = 5)
{
  if(class(submaps)[1] != "submapsList")
    stop("Need a submapsList matrix.")
  
  if(run_a_f) {
    submaps@likelihood_summary <- submapLikelihood(submaps@atlas)
    submaps@estimation_summary <- submapEstim(submaps@atlas)
    submaps@marker_summary <- summaryMarker(submaps@atlas, submaps@bedmatrix)
    submaps@submap_summary <- suppressWarnings(submapSummary(submaps@atlas))
  }
  
  test <- which(submaps@bedmatrix@ped$pheno == 2)
  if(run_a_f & probs & length(test) == 0 & missing(list.id) ) {
    cat("Cannot computes HBD, FLOD score and HFLOD score without individuals\n")
    cat("with pheno = 2. Use setSummary with a list.id argument.\n")
    return(submaps)
  }
  
  if(run_a_f && probs) {
    l1 <- setHBDprob(submaps, list.id=list.id, quality=quality)
    submaps <- l1[[1]]
    submaps <- setFLOD(submaps=submaps, condition=l1[[2]], q=q)
   
    # test class of first submap to check if recap is ok 
    if(!(class(submaps@atlas[[1]])[1] == "snpsMatrix" & recap.by.segments)) {
      l2 <- recap(submaps, recap.by.segments=recap.by.segments, list.id=list.id)
      submaps@HBD_recap <- l2[[1]]
      submaps@FLOD_recap <- l2[[2]]
      submaps@HBDsegments <- HBDsegments(submaps, threshold=threshold, n.consecutive.marker=n.consecutive.marker, recap.by.segments = recap.by.segments) 
      submaps@HFLOD <- setHFLOD(submaps)
    } else {
      warning("Cannot run HBD.recap, FLOD.recap and HFLOD with recap.by.segments = TRUE when submaps are created by 'Distance'")
    }
  }
  submaps
}


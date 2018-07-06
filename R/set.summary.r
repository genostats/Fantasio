#' Creation summary statistic files
#' 
#' This function is uses to ouput all the summary files after the creation of the submaps by the following two functions : 
#' `makeSubmapsByHotspots` and `makeSubmapsBySnps` is over.
#' 
#' 
#' @param submaps a list.submaps object
#' @param list.id you can either :
#'     - ignore this parameter if you want to compute HBD, FLOD and HFLOD 
#'       for individuals who are considerated INBRED and with a QUALITY
#'       greater or equal to 95%}
#'     - enter a list of individual for a computation of HBD, FLOD score HFLOD score for them
#'     - use "all" for a computation of HBD, FLOD score and HFLOD score for every individual
#' @param run_a_f a flag indicating if the estimation of a and f has been computed for the submaps using `festim` (default is TRUE)
#' @param probs a flag indicating if the HBD probabilities and FLOD score has been computed for the submaps (default is TRUE)
#' @param epsilon genotype error rate (default is 0.001)
#' @param run.festim whether you want to computes a, f, p.lrt, likelihood0/1 for each submaps  
#' @param run.proba whether you want to computes HBD, FLOD score and HFLOD score  
#' @param by_segments if the summary files has to be computed considering segments or snps (defaut is FALSE) 
#' (more information in the documentation of `makeSubmapsByHotspots` or`makeSubmapsBySnps` functions)
#' @param threshold the value of the threshold when finding HBD segment, threshold is the probability of being HBD or not (default is 0.5)
#' @param q  Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param quality The minimum percentage use to assume if a submap is valid (default is 95)
#' @param n.consecutive.marker the number of consecutive marker with a probabilitie equal or greater to the value of the threshold, to be use to fing HBD segments (default is 5)
#' @param test a vector containing the index of the different individual with a STATUS of 2 (individual with the disease)
#' 
#' @details the function filled the empty slots of the list.submaps object. It calls different functions and uses the results of each one to filled 
#' the object. This function computes summary for the following elements :
#' @details - summary for likelihood0/likelihood1
#' @details - summary for estimation of a and f
#' @details - marker summary
#' @details - summary for the submaps
#' @details - summary for HBD probabilities
#' @details - summary for HFLOD
#' @details - comutation of HBD segments
#' 
#' @return return a new list object containing every dataframe and object created 
#' 
#' @seealso \code{\link{makeSubmapsByHotspots}}
#' @seealso \code{\link{createSegmentsListByHotspots}}
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
#' submaps <- makeSubmapsByHotspots(bedMatrix, 5, segmentList) #this function is a wrapper that uses the function setSummary
#' 
#' @export
setSummary <- function(submaps, list.id, run_a_f = TRUE, probs = TRUE, by_segments=FALSE, q=1e-4, threshold=0.5, quality=95, n.consecutive.marker=5)
{
  if(class(submaps)[1] != "list.submaps")
    stop("Need a list.submaps matrix to eat")
   
  
  if(run_a_f)
  {
    submaps@likelihood_summary <- submapLikelihood(submaps@atlas)
  	submaps@estimation_summary <- submapEstim(submaps@atlas)
  	submaps@marker_summary <- summaryMarker(submaps@atlas, submaps@bedmatrix)
  	submaps@submap_summary <- submapSummary(submaps@atlas)
  }
  
  test <- which(submaps@bedmatrix@ped$pheno == 2)
  if(length(test) == 0)
  {
    cat("Cannot computes HBD, FLOD score and HFLOD score without 
        individuals with a STATUS 2, you can try using a the list.id arguments.\n")
    return(submaps)
  }
  
  if(run_a_f && probs)
  {
    l1 <- set.HBD.prob(submaps, list.id=list.id, quality=quality)
    submaps <- l1[[1]]
    submaps <- set.FLOD(submaps=submaps, condition=l1[[2]], q=q)
    
    if(!(class(submaps@atlas[[1]])[1] == "snps.matrix" & by_segments))
    {
      l2 <- recap(submaps, by_segments=by_segments, list.id=list.id)
  	  submaps@HBD_recap <- l2[[1]]
  	  submaps@FLOD_recap <- l2[[2]]
  	  submaps@HBD_segments <- HBD.segments(submaps, threshold=threshold, n.consecutive.marker=n.consecutive.marker, by_segments = by_segments) 
  	  submaps@HFLOD <- set.HFLOD(submaps)
    }else{
      cat("Cannot produce HBD.recap, FLOD.recap and HFLOD with option recap.by.segment to TRUE when using gap between markers (by distance),
      use hotposts with this option instead to make your submap\n")
    }
  	
  }  
  submaps
}
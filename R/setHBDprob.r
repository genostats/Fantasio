#' Computation of HBD probabilities
#' 
#' This function is used to compute HBD probabilities on individuals in a sample
#' 
#' @param submaps A atlas object
#' @param list.id you can either :
#'     - ignore this parameter if you want to compute HBD, FLOD and HFLOD 
#'       for individuals who are considerated INBRED and with a QUALITY
#'       greater or equal to 95%}
#'     - enter a list of individual for a computation of HBD, FLOD score HFLOD score for them
#'     - use "all" for a computation of HBD, FLOD score and HFLOD score for every individual
#' @param quality The minimum percentage use to assume if a submap is valid (default is 95)
#' @details This function iterates over the slots submaps_list of the atlas object.
#' @details For each submaps in the slots submaps_list of the object, the slot HBD.prob will be filled with a matrix of dimension : number_inidividual x number_of_markers
#' @details By default the function only computes HBD probabilities for INBRED individuals and with a quality equal or greater than 95%. However if you pass the keyword "all" to 
#' the list.id argument, this function will then computes HBD probabilities for all the individuals in your data. If you want a specific individual, then give a vector containing the family id 
#' and the individual id separated by an underscore to the list.id argument.
#' 
#' @return the atlas object with each HBD.prob slot of each submaps in the slot submaps_list computed
#' 
#' @seealso setFLOD
#' @seealso setHFLOD
#' 
#' @export
setHBDprob <- function(submaps, list.id, quality = 95)
{
  if(class(submaps@submaps_list[[1]])[1] != "snpsMatrix" & class(submaps@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  if(!missing(list.id))
  {
    if(list.id == "all")
    {
      condition <- 1:nrow(submaps@submap_summary)
    }else{
      vec <- strsplit(list.id, "_")
      condition <- sapply(vec, function(i) which(submaps@submap_summary$FID == i[1] & submaps@submap_summary$IID == i[2]))
    }
  }else{
    condition <- which(submaps@submap_summary$QUALITY >= quality & submaps@submap_summary$INBRED)
    if(length(condition) == 0) {
      warning(strwrap(paste("No inbred found with the following default parameters : QUALITY = ", quality , ", inbred = TRUE.
          You can try to change quality parameters or give a vector of individual with 'list.id' or use list.id = \"all\". 
          We will use all the individuals in the sample"), prefix = " ", initial=""))
      #for simplicity sake, not returning an empty results
      condition <- 1:nrow(submaps@submap_summary)
    }
      
  }
  
  id    <- as.vector(submaps@submap_summary$IID[condition])
  famid <- as.vector(submaps@submap_summary$FID[condition])
  
  for(i in 1:length(submaps@submaps_list))
  {
    HBD_prob <- matrix(NA, nrow = length(condition), ncol = submaps@submaps_list[[i]]@ncol)#HBD matrix
    dimnames(HBD_prob) <- list(rownames(HBD_prob) <- paste(famid,id, sep = "_"), colnames(HBD_prob) <- c(submaps@submaps_list[[i]]@map$id))
    
    for (j in 1:nrow(HBD_prob))
    {
      j1 <- condition[j]
      if(!is.na(submaps@submaps_list[[i]]@a[j1]) & (submaps@submaps_list[[i]]@a[j1] <= 1) & !is.na(submaps@submaps_list[[i]]@f[j1]))
      {
        HBD_prob[j,1:ncol(HBD_prob)] <-forward.backward(get.log.emiss(submaps@submaps_list[[i]], j1), 
                                                        submaps@submaps_list[[i]]@delta.dist, 
                                                        submaps@submaps_list[[i]]@a[j1],
                                                        submaps@submaps_list[[i]]@f[j1] )[2,]
      }
    }
    HBD_prob[!is.finite(HBD_prob)] <- 0
    submaps@submaps_list[[i]]@HBD.prob <- HBD_prob 
  }
  l <- list(submaps, condition)
  return(l)
}





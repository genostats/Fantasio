#' Estimation of likelihood0 and likelihood2
#' 
#' This function creates a summary on the submaps created. 
#' 
#' @param submaps a list of submaps 
#' 
#' @details the first element of the list is the estimation of likelihood0
#' @details the second element of the list is the estimation of likelihood1
#' 
#' @return this function returns a list of dataframe.
#' 
#' @seealso setSummary
#' 
#' 
#' @export

submapLikelihood <- function(submaps)
{
  if(class(submaps[[1]])[1] != "snps.matrix" & class(submaps[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.") 
  
  likelihood <- list()
  df0 <- data.frame(likelihoodH0 = sapply(submaps, function(x) x@likelihood0))
  df1 <- data.frame(likelihoodH1 = sapply(submaps, function(x) x@likelihood1))
  
  likelihood[[1]] <- df1
  likelihood[[2]] <- df0
  
  return(likelihood)       
}
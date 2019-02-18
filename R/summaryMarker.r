#' How many markers have been selected 0, 1, ... times
#' 
#' Gives the number of markers of the initial bedmatrix that have never been picked in any submap,
#' that have been picked once, etc.
#'
#' @param atlas an atlas 
#' 
#' @return A data frame with columns
#' \describe{
#'   \item{markers}{number of markers}
#'   \item{picked}{number of times picked}
#' } 
#' 
#' @seealso setSummary
#' 
#' @export
summaryMarker <- function(atlas) {
  submaps <- atlas@submaps_list
  bedmatrix <- atlas@bedmatrix
  
  b <- markerRepresentation(atlas)
  res <- numeric(length(submaps))

  for(i in 1:length(submaps)) {
    res[i] <- sum(b$Freq == i)
  }
  zero <- length(bedmatrix@snps$chr) - sum(res)
  
  df <- data.frame( markers = c(zero, res), picked = 0:length(res))
  df
}

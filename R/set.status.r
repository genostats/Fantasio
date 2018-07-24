#' Set status in the bed.matrix
#' 
#' Use this function if you want to change the status of your individuals in the bed.matrix
#' 
#' @param bedmatrix a bed.matrix 
#' @param status value accepted are 0, -9, 1, 2 (default is 2)
#' 
#' @export
set.status <- function(bedmatrix, status=2)
{
  if(status %in% c(0, -9, 1, 2))
    bedmatrix@ped$pheno <- rep(status, length(bedmatrix@ped$id))
  else 
    stop("values accepted for status are 0, -9, 1, 2")
  
  bedmatrix
}
#' Set p value in the bed.matrix
#' 
#' Use this function if you want to change the allelic frequencies in the bed.matrix object.
#' 
#' @param x a bed.matrix 
#' @param freq vector of allelic frequencies
#' 
#' @export
setAllelicFrequencies <- function(x, freq)
{
  if(length(freq) != ncol(x))
    stop("Length of 'freq' and number of SNPs in 'x' mismatch")
  
  x@p <- as.numeric(freq)
  x
}

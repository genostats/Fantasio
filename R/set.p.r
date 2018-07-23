#' Set p value in the bed.matrix
#' 
#' Use this function if you want to change the allelic frequencies in the bed.matrix object.
#' @param x a bed.matrix 
#' @param freq vector of allelic frequencies
#' 
#@export
set.allelic.frequencies <- function(x, freq)
{
  if(length(freq) != length(x@p))
    stop("length of the replacement vector for slot p mismatch")
  
  x@p <- as.numeric(freq)
}
#' Computation of HFLOD scores
#' 
#' This function is used to compute HFLOD scores on individuals in a sample
#' 
#' @param atlas a atlas object
#' @param w.id vector of indices to be considered in the FLOD_recap matrix
#' 
#' @details This function iterates over the slots submaps_list of the atlas object.
#' 
#' @return the atlas object with it's slot HFLOD filled with a matrix of dimensions : number_of_marker x 2. The first column is the value of HFLOD for the marker. 
#' The second value is the moving average probability.
#' 
#' @seealso setHBDProb
#' @seealso setFLOD
#'
#' @export
setHFLOD <- function(atlas, w.id)
{
  HFLOD <- getPosition(atlas, atlas@FLOD_recap)

  HFLOD_value <- numeric(nrow(HFLOD))
  ALPHA_value <- numeric(nrow(HFLOD))
  
  for (j in seq_len(nrow(HFLOD))) {
    # function h(alpha)
    h <- function(alpha)
      return(sum(log10( alpha*exp(atlas@FLOD_recap[w.id,j]*log(10))+(1-alpha) ), na.rm = TRUE))
    
    # optimisation of h(alpha) ; 
    res <- optimize( h, c(0,1), maximum = TRUE, tol = 0.001 )
    
    HFLOD_value[j] <- res$objective # HFLOD = h(alpha max)
    ALPHA_value[j] <- res$maximum   # alpha max 
  }
  HFLOD$HFLOD <- HFLOD_value
  HFLOD$ALPHA <- ALPHA_value

  atlas@HFLOD <- HFLOD
  atlas 
}

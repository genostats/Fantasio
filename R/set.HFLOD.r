#' Computation of HFLOD scores
#' 
#' This function is uses to compute HFLOD scores on individual present in a sample
#' 
#' @param submaps a list.submaps object
#' 
#' @details This function iterates over the slots atlas of the list.submaps object.
#' 
#' @return the list.submaps object with it's slot HFLOD filled with a matrix of dimensions : number_of_marker x 2. The first column is the value of HFLOD for the marker. 
#' The second value is the moving average probability.
#' 
#' @seealso \code{\link{set.HBD.prob}}
#' @seealso \code{\link{set.FLOD}}
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListByHotspots(bedMatrix)
#' individualList <- c("familyID0_individualID0", "familyID1_individualID2"), "familyID2_individualID2")
#' makeSubmapsByHotspots(bedMatrix, 10, segmentList, list.id=individualList)  #the function set.HFLOD is use inside the function setSummary of this function
#' @export
set.HFLOD <- function(submaps)
{
  HFLOD <- matrix(0.0, nrow = ncol(submaps@FLOD_recap), ncol = 2)
  
  for (j in 1:nrow(HFLOD))
  {
    # function h(alpha)
    h <- function(alpha)
      return(sum(log10( alpha*exp(submaps@FLOD_recap[,j]*log(10))+(1-alpha) ), na.rm = TRUE))

    #maximisation of h(alpha) ; 
    res <- optimize( h, c(0,1), maximum = TRUE, tol = 0.001 )
    
    
    HFLOD[j,1]<-res$objective # HFLOD = h(alpha max)
    HFLOD[j,2]<-res$maximum  #alpha max 
  }
  return(HFLOD)
}


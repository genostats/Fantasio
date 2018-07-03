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
  if(class(submaps@atlas[[1]])[1] != "snps.matrix" & class(submaps@atlas[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.") 
  
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
  {
    stop("Need a bed.matrix to eat")
  }

  if(submaps@bySegments)
  {
    index <- sapply(submaps@atlas, function(i) i@submap) #index of the marker
    
    poscM <- matrix(0.0, nrow = nrow(index), ncol = ncol(index))
    posBp <- matrix(0.0, nrow = nrow(index), ncol = ncol(index))
    
    poscM_mean <- c()
    posBp_mean <- c()
    
    for(i in 1:nrow(index))
    {
      for(j in 1:ncol(index))#get position of the marker
      {
        poscM[i, j] <- x@snps$dist[index[i,j]]
        posBp[i, j] <- x@snps$pos[index[i,j]]
      }
      #calculate mean value
      poscM_mean <- c(poscM_mean, mean(poscM[i,]))
      posBp_mean <- c(posBp_mean, mean(posBp[i,]))
    }
    
    
    HFLOD <- data.frame(CHR = submaps@atlas[[1]]@map$chr, 
                        pos_cM = poscM_mean, 
                        pos_Bp = posBp_mean)
    
  }else{
    names <- colnames(submaps@FLOD_recap)
    index <- match(names, submaps@bedmatrix@snps$id)
    distance_cM <- submaps@bedmatrix@snps$dist[index]
    distance_bP <- submaps@bedmatrix@snps$pos[index]
    chr <- x@snps$chr[index]
    
    HFLOD <- data.frame(CHR = chr, 
                        SNPS = names, 
                        pos_cM = distance_cM,
                        pos_Bp =distance_bP)
  }
  # HFLOD <- matrix(0.0, nrow = ncol(submaps@FLOD_recap), ncol = 2)
  HFLOD_value <- c()
  ALPHA_value <- c()
  for (j in 1:nrow(HFLOD))
  {
    # function h(alpha)
    h <- function(alpha)
      return(sum(log10( alpha*exp(submaps@FLOD_recap[,j]*log(10))+(1-alpha) ), na.rm = TRUE))

    # optimisation of h(alpha) ; 
    res <- optimize( h, c(0,1), maximum = TRUE, tol = 0.001 )
    
    
    HFLOD_value <- c(HFLOD_value, res$objective) # HFLOD = h(alpha max)
    ALPHA_value <- c(ALPHA_value, res$maximum)   # alpha max 
  }
  HFLOD$HFLOD <- HFLOD_value
  HFLOD$ALPHA <- ALPHA_value
  return(HFLOD)
}

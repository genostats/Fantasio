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
#' @seealso setHBDprob
#' @seealso setFLOD
#'
#' @export
setHFLOD <- function(atlas, w.id)
{
  if(class(atlas@submaps_list[[1]])[1] != "snpsMatrix" & class(atlas@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  if(class(atlas@bedmatrix)[1] != "bed.matrix") 
    stop("Need a bed.matrix.")
  
  if(atlas@bySegments) {
    index <- sapply(atlas@submaps_list, function(i) i@submap) #index of the marker
    
    poscM <- matrix(0.0, nrow = nrow(index), ncol = ncol(index))
    posBp <- matrix(0.0, nrow = nrow(index), ncol = ncol(index))
    
    poscM_mean <- numeric(nrow(index))
    posBp_mean <- numeric(nrow(index))
    
    for(i in 1:nrow(index)) {
      for(j in 1:ncol(index))#get position of the marker
      {
        poscM[i, j] <- atlas@bedmatrix@snps$dist[index[i,j]]
        posBp[i, j] <- atlas@bedmatrix@snps$pos[index[i,j]]
      }
      #calculate mean value
      poscM_mean[i] <- mean(poscM[i,])
      posBp_mean[i] <- mean(posBp[i,])
    }
    
    HFLOD <- data.frame(CHR    = atlas@submaps_list[[1]]@map$chr, 
                        pos_cM = poscM_mean, 
                        pos_Bp = posBp_mean)
    
  } else {
    names       <- colnames(atlas@FLOD_recap)  # les ids des marqueurs où on a le FLOD
    index       <- match(names, atlas@bedmatrix@snps$id)
    distance_cM <- atlas@bedmatrix@snps$dist[index]
    distance_bP <- atlas@bedmatrix@snps$pos[index]
    chr         <- atlas@bedmatrix@snps$chr[index]
    
    HFLOD <- data.frame(CHR    = chr, 
                        SNPS   = names, 
                        pos_cM = distance_cM,
                        pos_Bp = distance_bP)
  }

  HFLOD_value <- numeric(nrow(HFLOD))
  ALPHA_value <- numeric(nrow(HFLOD))
  
  for (j in 1:nrow(HFLOD)) {
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

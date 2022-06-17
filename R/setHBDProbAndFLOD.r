#' Computation of HBD probabilities and FLOD scores by segments 
#' 
#' This function is used to compute HBD probabilities and FLOD scores on individuals in a sample
#' 
#' @param atlas an atlas object
#' @param w.id The indices of individals for which HBD probabilities are computed
#' @details This function iterates over the slots submaps_list of the atlas object.
#' @details For each submaps in the slots submaps_list of the object, the slots HBD_recap and FLOD_recap will be updated with a matrix of dimension : number_individual x number_of_segments
#' 
#' @return the atlas object with HBD_recap and FLOD_recap by segments
#' 
#' @seealso \code{\link{setHBDProbAndFLODBySnps}}
#' @seealso \code{\link{setHFLOD}}
#' 
#' @export

setHBDProbAndFLOD <- function(atlas, w.id, q = 1e-4)
{
  if(class(atlas@submaps_list[[1]])[1] != "snpsMatrix" & class(atlas@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  if(class(atlas@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  id    <- as.vector(atlas@submap_summary$id[w.id])
  famid <- as.vector(atlas@submap_summary$famid[w.id])
  
moy_HBD_prob <- matrix(0.0, nrow = length(w.id), ncol = atlas@submaps_list[[1]]@ncol)#HBD matrix
dimnames(moy_HBD_prob) <- list( uniqueIds(famid,id), lapply(c(1:atlas@submaps_list[[1]]@ncol), function(i) paste0('s', i)) )

moy_FLOD <- matrix(0.0, nrow = length(w.id), ncol = atlas@submaps_list[[1]]@ncol)#HBD matrix
dimnames(moy_FLOD) <- list( uniqueIds(famid,id), lapply(c(1:atlas@submaps_list[[1]]@ncol), function(i) paste0('s', i)) )  
  
for(i in seq_along(atlas@submaps_list)) {
  HBD_prob <- matrix(NA, nrow = length(w.id), ncol = atlas@submaps_list[[i]]@ncol)#HBD matrix
  dimnames(HBD_prob) <- list( uniqueIds(famid,id), atlas@submaps_list[[i]]@map$id)
  
  FLOD <- matrix(NA, nrow = length(w.id), ncol = atlas@submaps_list[[i]]@ncol)#HBD matrix
  dimnames(FLOD) <- list( uniqueIds(famid,id), atlas@submaps_list[[i]]@map$id)
  
  for (j in seq_len(nrow(HBD_prob))) {
    j1 <- w.id[j]
    if(!is.na(atlas@submaps_list[[i]]@a[j1]) & (atlas@submaps_list[[i]]@a[j1] <= 1) & !is.na(atlas@submaps_list[[i]]@f[j1])) {
      HBD_prob[j,seq_len(ncol(HBD_prob))] <-forwardBackward(getLogEmiss(atlas@submaps_list[[i]], j1), 
                                                            atlas@submaps_list[[i]]@delta.dist, 
                                                            atlas@submaps_list[[i]]@a[j1],
                                                            atlas@submaps_list[[i]]@f[j1] )[2,]
    }
    HBD_prob[!is.finite(HBD_prob)] <- 0
    
    if( (atlas@submaps_list[[i]]@a[w.id[j]] < 1) & !is.na(atlas@submaps_list[[i]]@f[w.id[j]]) ) {
      FLOD[j,] <- log10((HBD_prob[j,] + q * ( 1 - HBD_prob[j,]))/
                          (atlas@submaps_list[[i]]@f[w.id[j]] + q * ( 1 - atlas@submaps_list[[i]]@f[w.id[j]]))) 
    }
  }
  moy_HBD_prob <- moy_HBD_prob + HBD_prob/length(atlas@submaps_list)
  moy_FLOD <- moy_FLOD + FLOD/length(atlas@submaps_list)
}
atlas@HBD_recap <- moy_HBD_prob
atlas@FLOD_recap <- moy_FLOD
atlas
}

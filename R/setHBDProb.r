#' Computation of HBD probabilities
#' 
#' This function is used to compute HBD probabilities on individuals in a sample
#' 
#' @param atlas A atlas object
#' @param w.id The indices of individals for which HBD probabilities are computed
#' @details This function iterates over the slots submaps_list of the atlas object.
#' @details For each submaps in the slots submaps_list of the object, the slot HBD.prob will be filled with a matrix of dimension : number_individual x number_of_markers
#' 
#' @return the atlas object with each HBD.prob slot of each submaps in the slot submaps_list computed
#' 
#' @seealso setFLOD
#' @seealso setHFLOD
#' 
#' @export
setHBDProb <- function(atlas, w.id)
{
  if(class(atlas@submaps_list[[1]])[1] != "snpsMatrix" & class(atlas@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  if(class(atlas@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  id    <- as.vector(atlas@submap_summary$id[w.id])
  famid <- as.vector(atlas@submap_summary$famid[w.id])
  
  for(i in seq_along(atlas@submaps_list)) {
    HBD_prob <- matrix(NA, nrow = length(w.id), ncol = atlas@submaps_list[[i]]@ncol)#HBD matrix
    dimnames(HBD_prob) <- list( UniqueIds(famid,id), atlas@submaps_list[[i]]@map$id)
    
    for (j in seq_len(nrow(HBD_prob))) {
      j1 <- w.id[j]
      if(!is.na(atlas@submaps_list[[i]]@a[j1]) & (atlas@submaps_list[[i]]@a[j1] <= 1) & !is.na(atlas@submaps_list[[i]]@f[j1])) {
        HBD_prob[j,seq_len(ncol(HBD_prob))] <-forward.backward(get.log.emiss(atlas@submaps_list[[i]], j1), 
                                                               atlas@submaps_list[[i]]@delta.dist, 
                                                               atlas@submaps_list[[i]]@a[j1],
                                                               atlas@submaps_list[[i]]@f[j1] )[2,]
      }
    }
    HBD_prob[!is.finite(HBD_prob)] <- 0
    atlas@submaps_list[[i]]@HBD.prob <- HBD_prob 
  }
  atlas
}




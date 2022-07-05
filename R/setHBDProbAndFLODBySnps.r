#' Computation of HBD probabilities and FLOD scores by snps 
#' 
#' This function is used to compute HBD probabilities and FLOD scores on individuals in a sample
#' 
#' @param atlas an atlas object
#' @param w.id The indices of individals for which HBD probabilities are computed
#' @details This function iterates over the slots submaps_list of the atlas object.
#' @details For each submaps in the slots submaps_list of the object, the slots HBD_recap and FLOD_recap will be updated with a matrix of dimension : number_individual x number_of_markers
#' 
#' @return the atlas object with HBD_recap and FLOD_recap by snps
#' 
#' @seealso \code{\link{setHFLOD}}
#' 
#' @export

setHBDProbAndFLODBySnps <- function(atlas, w.id, q = 1e-4)
{
  if(class(atlas@submaps_list[[1]])[1] != "snpsMatrix" & class(atlas@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  if(class(atlas@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  id    <- as.vector(atlas@submap_summary$id[w.id])
  famid <- as.vector(atlas@submap_summary$famid[w.id])
  
  h <- new.env()
  
  for(i in seq_along(atlas@submaps_list)) {
    HBD_prob <- matrix(0.0, nrow = length(w.id), ncol = atlas@submaps_list[[i]]@ncol)#HBD matrix
    dimnames(HBD_prob) <- list( uniqueIds(famid,id), atlas@submaps_list[[i]]@map$id)
    
    FLOD <- matrix(0.0, nrow = length(w.id), ncol = atlas@submaps_list[[i]]@ncol)#HBD matrix
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
    
      if( (atlas@submaps_list[[i]]@a[j1] < 1) & !is.na(atlas@submaps_list[[i]]@f[j1]) ) {
        FLOD[j,] <- log10((HBD_prob[j,] + q * ( 1 - HBD_prob[j,]))/
                          (atlas@submaps_list[[i]]@f[j1] + q * ( 1 - atlas@submaps_list[[i]]@f[j1]))) 
      }
    }

    for (k in seq_along(atlas@submaps_list[[i]]@map$id)) {
      if (length(h[[atlas@submaps_list[[i]]@map$id[[k]]]][[1]]) == 0) {
        h[[atlas@submaps_list[[i]]@map$id[[k]]]][['HBD_prob']] <- HBD_prob[,k]
        h[[atlas@submaps_list[[i]]@map$id[[k]]]][['FLOD']] <- FLOD[,k]
        h[[atlas@submaps_list[[i]]@map$id[[k]]]][['counts']] <- 1
      } else { 
        h[[atlas@submaps_list[[i]]@map$id[[k]]]][['HBD_prob']] <- h[[atlas@submaps_list[[i]]@map$id[[k]]]][['HBD_prob']] + HBD_prob[,k]
        h[[atlas@submaps_list[[i]]@map$id[[k]]]][['FLOD']] <- h[[atlas@submaps_list[[i]]@map$id[[k]]]][['FLOD']] + FLOD[,k]
        h[[atlas@submaps_list[[i]]@map$id[[k]]]][['counts']] <- h[[atlas@submaps_list[[i]]@map$id[[k]]]][['counts']] + 1
      }
    }
  }
  
  for (n in names(h)) {
  h[[n]][['HBD_prob']] <- h[[n]][['HBD_prob']] / h[[n]][['counts']]
  h[[n]][['FLOD']] <- h[[n]][['FLOD']] / h[[n]][['counts']]
  }
  
  hbd <- NULL
  hbd <- as.data.frame(lapply(names(h), function (i) cbind(hbd, h[[i]][['HBD_prob']])))
  colnames(hbd) <- names(h)
  hbd <- as.matrix(hbd)
  
  atlas@HBD_recap <- hbd
  
  flod <- NULL
  flod <- as.data.frame(lapply(names(h), function (i) cbind(flod, h[[i]][['FLOD']])))
  colnames(flod) <- names(h)
  flod <- as.matrix(flod)
  
  atlas@FLOD_recap <- flod
  
  po <- match(colnames(atlas@HBD_recap), atlas@bedmatrix@snps$id)
  snp.chr <- atlas@bedmatrix@snps$chr[po]
  snp.pos <- atlas@bedmatrix@snps$pos[po]
  atlas@HBD_recap <- atlas@HBD_recap[, order(snp.chr, snp.pos)]
  atlas@FLOD_recap <- atlas@FLOD_recap[, order(snp.chr, snp.pos)]
  
  atlas
}

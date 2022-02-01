getPosition = function (x, hbd) {

	if(class(x@submaps_list[[1]])[1] != "snpsMatrix" & class(x@submaps_list[[1]])[1] != "HostspotsMatrix")
    		stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  	if(class(x@bedmatrix)[1] != "bed.matrix") 
    		stop("Need a bed.matrix.")
  
  	# create the structure of the final dataframe
  	# recovery of chr, snps, pos_cM, pos_Bp
  	if(x@bySegments) {
    		index <- sapply(x@submaps_list, function(i) i@submap) #index of the marker
    
    		poscM <- matrix(0.0, nrow = nrow(index), ncol = ncol(index))
    		posBp <- matrix(0.0, nrow = nrow(index), ncol = ncol(index))
    
    		poscM_mean <- numeric(nrow(index))
    		posBp_mean <- numeric(nrow(index))
    
    	for(i in seq_len(nrow(index))) {
      		for(j in seq_len(ncol(index)))#get position of the marker
      		{
       		poscM[i, j] <- x@bedmatrix@snps$dist[index[i,j]]
        		posBp[i, j] <- x@bedmatrix@snps$pos[index[i,j]]
      		}
      		#calculate mean value
      		poscM_mean[i] <- mean(poscM[i,])
      		posBp_mean[i] <- mean(posBp[i,])
   	 }
    
    	final <- data.frame(chr    = x@submaps_list[[1]]@map$chr, 
                            pos_cM = poscM_mean, 
                            pos_Bp = posBp_mean)
  	} else {
    		names       <- colnames(hbd)  # les ids des marqueurs où on a le FLOD/ la proba HBD
    		index       <- match(names, x@bedmatrix@snps$id)
    		distance_cM <- x@bedmatrix@snps$dist[index]
    		distance_bP <- x@bedmatrix@snps$pos[index]
    		chr         <- x@bedmatrix@snps$chr[index]
    
    	final <- data.frame(chr    = chr, 
                            snps   = names, 
                            pos_cM = distance_cM,
                            pos_Bp = distance_bP)
	}
}

##################################################################################
#This function creates HBD and FLOD recap dataframe by snps                      #
#                                                                                #
#!!! submaps : the list of object                                                #                                       
#!!! proba_HBD : a list of dataframe with the HBD probabilities                  #
#!!! proba_FLOD : a list of FLOD score                                           #
#                                                                                # 
#*** return a list of two dataframes with HBD and FLOD recap by snps in it       #
##################################################################################

recap.by.snps <- function(submaps, proba_HBD, proba_FLOD) {
  ##count the number of times the marker has been picked
  marqueurs <- as.data.frame(table(unlist(sapply(proba_HBD, function(x) colnames(x)))), stringsAsFactors=FALSE)
  
  ##count the number of times an individuals has appeared
  individuals <- as.vector(rownames(submaps@submaps_list[[1]]@HBD.prob))
  

  matrice_HBD <- matrix(0, nrow = length(individuals), ncol = nrow(marqueurs))
  matrice_FLOD <- matrix(0, nrow = length(individuals), ncol = nrow(marqueurs))
  
  nb_probs <- matrix(0L, nrow = length(individuals), ncol = nrow(marqueurs))
  
  dimnames(matrice_HBD)  <- list(individuals, marqueurs[,1])
  
  #to order chromosome in the matrix
 
  
  po <- match(colnames(matrice_HBD), submaps@bedmatrix@snps$id)
  snp.chr <- submaps@bedmatrix@snps$chr[po]
  snp.pos <- submaps@bedmatrix@snps$pos[po]
  matrice_HBD <- matrice_HBD[, order(snp.chr, snp.pos)]
  
  dimnames(matrice_FLOD) <- list(individuals, colnames(matrice_HBD))
    
  #find the index of the marker in every submaps
  az <- vapply(proba_HBD, function(jj) match(colnames(matrice_HBD), colnames(jj)), integer(length(marqueurs[,1])))

  
  for( i in 1:ncol(matrice_HBD)) {
    az1 <- az[i,]
    I <- which(!is.na(az1))
    v <- sapply(I, function(i) proba_HBD[[i]][, az1[i]], simplify = F)
    w <- sapply(I, function(i) proba_FLOD[[i]][, az1[i]], simplify = F)
    for( j in 1:length(v)) {
      matrice_HBD[,i] <- matrice_HBD[,i] + v[[j]]
      matrice_FLOD[,i] <- matrice_FLOD[,i] + w[[j]]
      nb_probs[,i] <- nb_probs[,i] + 1
    }
  }
  matrice_HBD  <- matrice_HBD/nb_probs
  matrice_FLOD <- matrice_FLOD/nb_probs
  l <- list(matrice_HBD, matrice_FLOD)
  return(l)
}

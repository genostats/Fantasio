##################################################################################
#This function creates HBD and FLOD recap dataframe                              #
#                                                                                #
#!!! submaps : the list of object                                                #                                       
#!!! proba_HBD : a list of dataframe with the HBD probabilities                  #
#!!! FLOD : a list of dataframe with the FLOD score                              #
#                                                                                # 
#*** return a list of two dataframes with HBD and FLOD recap in it               #
##################################################################################


recap.by.segments <- function(submaps, proba_HBD, FLOD) {

  marker_names <- colnames(submaps@submaps_list[[1]]@HBD.prob) # markers names
  correspondance <- match(marker_names, submaps@bedmatrix@snps$id)
  
  chr <- submaps@bedmatrix@snps$chr[correspondance] # chromosome corresponding to this marker
  
  ##trop lourd ?
  columns_names <- paste(rep("Segment",length(marker_names)), seq(1,length(marker_names)), rep("chr", length(marker_names)), chr, sep = "_")
  
  
  nom    <- as.vector(rownames(submaps@submaps_list[[1]]@HBD.prob)) # get the name of individuals
  matrice_HBD <- matrix(0, nrow=length(nom), ncol=ncol(proba_HBD[[1]]))
  matrice_FLOD <- matrix(0, nrow=length(nom), ncol=ncol(FLOD[[1]]))
  
  rownames(matrice_HBD) <- nom
  colnames(matrice_HBD) <- columns_names
  
  rownames(matrice_FLOD) <- nom
  colnames(matrice_FLOD) <- columns_names
  
  # loop over the individuals
  for(j in 1 : length(nom)) {
    Sum_HBD <- 0
    Sum_FLOD <- 0
    cpt <- 0
    
    for(i in seq_along(proba_HBD))   # every submaps
    {
      # get the line corresponding to our individual in the submap i 
      line <- which(nom[j] == rownames(proba_HBD[[i]]))  
      if(length(line) == 0) next()
      Sum_HBD <- Sum_HBD + proba_HBD[[i]][line, ] # sum every line
      
      line <- which(nom[j] == rownames(FLOD[[i]]))
      if(length(line) == 0) next()
      Sum_FLOD <- Sum_FLOD + FLOD[[i]][line, ]
      
      cpt <- cpt + 1 
    }
    matrice_HBD[j, ]  <- Sum_HBD  / cpt   # mean
    matrice_FLOD[j, ] <- Sum_FLOD / cpt
  }
  l <- list(matrice_HBD, matrice_FLOD)
  return(l)
}


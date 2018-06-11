##################################################################################
#This function creates HBD and FLOD recap dataframe                              #
#                                                                                #
#!!! submaps : the list of object                                                #                                       
#!!! proba_HBD : a list of dataframe with the HBD probabilities                  #
#!!! proba_FLOD : a list of dataframe with the FLOD score                        #
#                                                                                # 
#*** return a list of two dataframes with HBD and FLOD recap in it               #
##################################################################################


recap.by.segments <- function(submaps, proba_HBD, proba_FLOD)
{

  #max <- sapply(submaps@atlas, function(h) h@ncol)
  #max <- which(max == max(max))
  
  marker_names <- colnames(submaps@atlas[[1]]@HBD.prob)#recuperer le nom des marqueurs
  correspondance <- match(marker_names, submaps@bedmatrix@snps$id)
  
  chr <- submaps@bedmatrix@snps$chr[correspondance]#a quel chromosome appartient ce marqueur
  
  ##trop lourd ?
  columns_names <- paste(rep("Segment",length(marker_names)), seq(1,length(marker_names)), rep("chr", length(marker_names)), chr, sep = "_")
  
  
  nom    <- as.vector(rownames(submaps@atlas[[1]]@HBD.prob)) #recuperer le nom de chaque individus apparut
  #family <- as.vector(submaps@submap_summary$FID[which(submaps@submap_summary$INBRED == TRUE)])
  
  matrice_HBD <- matrix(0, nrow=length(nom), ncol=ncol(proba_HBD[[1]]))#creer une matrice avec les individus en lignes 
  matrice_FLOD <- matrix(0, nrow=length(nom), ncol=ncol(proba_FLOD[[1]]))
  
  #et les segments en colonnes
  rownames(matrice_HBD) <- nom
  colnames(matrice_HBD) <- columns_names
  
  rownames(matrice_FLOD) <- nom
  colnames(matrice_FLOD) <- columns_names
  
  #boucle de remplissage
  for(j in 1 : length(nom))#parcourir chaque individus
  {
    Sum_HBD <- 0
    Sum_FLOD <- 0
    cpt <- 0
    
    for(i in 1:length(proba_HBD))#parcourir chaque sous-cartes 
    {
      line <- which(nom[j] == rownames(proba_HBD[[i]]))#recuperer la lignes correspondant a notre individu pour la sous-carte i
      if(length(line) == 0) next()#si l'individu n'est pas dans la sous-carte on passe
      Sum_HBD <- Sum_HBD + proba_HBD[[i]][line, ]#on somme chaque ligne entres elles -> peut y avoir different marqueur mais meme segment
      
      line <- which(nom[j] == rownames(proba_FLOD[[i]]))#recuperer la lignes correspondant a notre individu pour la sous-carte i
      if(length(line) == 0) next()#si l'individu n'est pas dans la sous-carte on passe
      Sum_FLOD <- Sum_FLOD + proba_FLOD[[i]][line, ]#on somme chaque ligne entres elles -> peut y avoir different marqueur mais meme segment
      
      cpt <- cpt + 1 
      
    }
    matrice_HBD[j, ]  <- Sum_HBD  / cpt#on fait la moyenne
    matrice_FLOD[j, ] <- Sum_FLOD / cpt
  }
  l <- list(matrice_HBD, matrice_FLOD)
  return(l)
}


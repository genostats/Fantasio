##################################################################################
#This function creates HBD and FLOD recap dataframe by snps                      #
#                                                                                #
#!!! submaps : the list of object                                                #                                       
#!!! proba_HBD : a list of dataframe with the HBD probabilities                  #
#!!! proba_FLOD : a list of FLOD score                                           #
#                                                                                # 
#*** return a list of two dataframes with HBD and FLOD recap by snps in it       #
##################################################################################

recap.by.snps <- function(submaps, proba_HBD, proba_FLOD)
{
  ##tableau comptant le nombre de fois ou un marqueur a ete selectionne
  marqueurs <- as.data.frame(table(unlist(sapply(proba_HBD, function(x) colnames(x)))), stringsAsFactors=FALSE)
  
  ##tableau comptant le nombre de fois ou un individus est apparu
  individuals <- as.vector(rownames(submaps@atlas[[1]]@HBD.prob))
  
  #creation de la matrice = individus x marqueurs
  matrice_HBD <- matrix(0, nrow = length(individuals), ncol = nrow(marqueurs))
  matrice_FLOD <- matrix(0, nrow = length(individuals), ncol = nrow(marqueurs))
  
  nb_probs <- matrix(0L, nrow = length(individuals), ncol = nrow(marqueurs))
  
  dimnames(matrice_HBD)  <- list(individuals, marqueurs[,1])
  
  #to order chromosome in the matrix
 
  
  po <- match(colnames(matrice_HBD), submaps@bedmatrix@snps$id)
  snp.chr <- submaps@bedmatrix@snps$chr[po]
  snp.pos <- submaps@bedmatrix@snps$pos[po]
  #voo <- order(vi)
  matrice_HBD <- matrice_HBD[, order(snp.chr, snp.pos)]
  
  dimnames(matrice_FLOD) <- list(individuals, colnames(matrice_HBD))
    
  
  #trouver l'/les indice(s) et la/les sous-carte(s) avec des valeurs pour chaque marqueur
  az <- vapply(proba_HBD, function(jj) match(colnames(matrice_HBD), colnames(jj)), integer(length(marqueurs[,1])))
  
  #nom <- strsplit(rownames(matrice_HBD), "_")
  #nom <- sapply(nom, function(name) name[1])
  #M <- sapply(proba_HBD, function(pr) match(rownames(pr), rownames(matrice_HBD)))#trouver l'indice des individus apparuent dans chaque sous-cartes
  
  
  for( i in 1:ncol(matrice_HBD))
  {
    az1 <- az[i,]
    I <- which(!is.na(az1))
    v <- sapply(I, function(i) proba_HBD[[i]][, az1[i]], simplify = F)
    w <- sapply(I, function(i) proba_FLOD[[i]][, az1[i]], simplify = F)
    for( j in 1:length(v))
    {
      #m <- M[[ I[j] ]]
      #right_index <- which(!is.na(m))
      #m <- m[!is.na(m)]
      #matrice_HBD[m,i] <- matrice_HBD[m,i] + v[[j]][right_index]
      #matrice_FLOD[m,i] <- matrice_FLOD[m,i] + w[[j]][right_index]
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
HBD.recap <- function(submaps, by_segments=F)
{
  proba <- submap.HBD(submaps)#recuperer les probas HBD pour chaque sous cartes sous forme d'une liste

  if(by_segments)
  {
    #donner un nom aux colonnes
    marker_names <- colnames(submaps@atlas[[1]]@HBD.prob)#recuperer le nom des marqueurs
    correspondance <- match(marker_names, submaps@bedmatrix@snps$id)
    
    chr <- submaps@bedmatrix@snps$chr[correspondance]#a quel chromosome appartient ce marqueur
    
    ##trop lourd ?
    columns_names <- paste(rep("Segment",length(marker_names)), seq(1,length(marker_names)), rep("chr", length(marker_names)), chr, sep = "_")
    
    
    nom <- as.vector(submaps@submap_summary[which(submaps@submap_summary$INBRED == TRUE),]$IID) #recuperer le nom de chaque individus apparut
    
    matrice <- matrix(NA, nrow=length(nom), ncol=ncol(proba[[1]]))#creer une matrice avec les individus en lignes 
                                                                  #et les segments en colonnes
    rownames(matrice) <- nom
    colnames(matrice) <- columns_names
    
    #boucle de remplissage
    for(j in 1 : length(nom))#parcourir chaque individus
    {
      Sum <- 0
      cpt <- 0
      for(i in 1:length(proba))#parcourir chaque sous-cartes 
      {
        line <- which(nom[j] == rownames(proba[[i]]))#recuperer la lignes correspondant a notre individu pour la sous-carte i
        if(length(line) == 0) next()#si l'individu n'est pas dans la sous-carte on passe
        Sum <- Sum + proba[[i]][line, ]#on somme chaque ligne entres elles -> peut y avoir different marqueur mais meme segment
        cpt <- cpt + 1 
      }
      matrice[j, ] <- Sum / cpt#on fait la moyenne
      
    }
    return(matrice)
  }
  
  ##tableau comptant le nombre de fois ou un marqueur a ete selectionne
  marqueurs <- as.data.frame(table(unlist(sapply(proba, function(x) colnames(x)))), stringsAsFactors=FALSE)
  
  ##tableau comptant le nombre de fois ou un individus est apparu
  individuals <- as.vector(submaps@submap_summary[which(submaps@submap_summary$INBRED == TRUE),]$IID)
  
  #creation de la matrice = individus x marqueurs
  matrice <- matrix(0, nrow = length(individuals), ncol = nrow(marqueurs))
  nb_probs <- matrix(0L, nrow = length(individuals), ncol = nrow(marqueurs))
  
  dimnames(matrice) <- list(individuals, marqueurs[,1])
  
  #to order chromosome in the matrix
  
  po <- match(colnames(matrice), submaps@bedmatrix@snps$id)
  snp.chr <- submaps@bedmatrix@snps$chr[po]
  snp.pos <- submaps@bedmatrix@snps$pos[po]
  #voo <- order(vi)
  matrice <- matrice[, order(snp.chr, snp.pos)]
  
  
  #trouver l'/les indice(s) et la/les sous-carte(s) avec des valeurs pour chaque marqueur
  az <- vapply(proba, function(jj) match(marqueurs[,1], colnames(jj)), integer(length(marqueurs[,1])))
  
 
  M <- sapply(proba, function(pr) match(rownames(pr), rownames(matrice)))
  
  
  for( i in 1:ncol(matrice))
  {
    az1 <- az[i,]
    I <- which(!is.na(az1))
    v <- sapply(I, function(i) proba[[i]][, az1[i]], simplify = F)
    for( j in 1:length(v))
    {
      m <- M[[ I[j] ]]
      right_index <- which(!is.na(m))
      m <- m[!is.na(m)]
      matrice[m,i] <- matrice[m,i] + v[[j]][right_index]
      nb_probs[m,i] <- nb_probs[m,i] + 1
    }
  }
  matrice/nb_probs
}  

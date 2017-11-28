HBD.recap <- function(submaps, by_segments=F)
{
  proba <- submap.HBD(submaps)#recuperer les probas HBD pour chaque sous cartes sous forme d'une liste

  if(by_segments)
  {
    #donner un nom aux colonnes
    marker_names <- colnames(submaps@atlas[[1]]@HBD.prob)#recuperer le nom des marqueurs
    chr <- submaps@bedmatrix@snps$chr[match(marker_names, submaps@bedmatrix@snps$id)]#a quel chromosome appartient ce marqueur
    columns_names <- paste(rep("Segment",length(marker_names)), seq(1,length(marker_names)), rep("chr", length(marker_names)), chr, sep = "_")
    
    
    nom <- unique(unlist(sapply(proba, function(x) rownames(x)))) #recuperer le nom de chaque individus apparut
    
    matrice <- matrix(NA, nrow=length(nom), ncol=ncol(proba[[1]]))#creer une matrice avec les individus en lignes 
                                                                  #et les segments en colonnes
    rownames(matrice) <- nom
    colnames(matrice) <- columns_names
    
    #boucle de remplissage
    Sum <- 0
    for(j in 1 : length(nom))#parcourir chaque individus
    {
      cpt <- 0
      for(i in 1:length(proba))#parcourir chaque sous-cartes 
      {
        line <- which(nom[j] == rownames(proba[[i]]))#recuperer la lignes correspondant a notre individu
        if(length(line) == 0) next()#si l'individu n'est pas dans la sous-carte on passe
        Sum <- Sum + proba[[i]][line, ]#on somme chaque ligne entres elles
        cpt <- cpt + 1 
      }
      matrice[j, ] <- Sum / cpt#on fait la moyenne
      
    }
    return(matrice)
  }
  
  ##tableau comptant le nombre de fois ou un marqueur a ete selectionne
  marqueurs <- as.data.frame(table(unlist(sapply(proba, function(x) colnames(x)))), stringsAsFactors=FALSE)
  
  ##tableau comptant le nombre de fois ou un individus est apparu
  individuals <- as.data.frame(table(unlist(sapply(proba, function(x) rownames(x)))), stringsAsFactors=FALSE) 

  
  #creation de la matrice = individus x marqueurs
  matrice <- matrix(0, nrow = nrow(individuals), ncol = nrow(marqueurs))
  nb_probs <- matrix(0L, nrow = nrow(individuals), ncol = nrow(marqueurs))
  
  dimnames(matrice) <- list(individuals[,1], marqueurs[,1])
  
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
      matrice[m,i] <- matrice[m,i] + v[[j]]
      nb_probs[m,i] <- nb_probs[m,i] + 1
    }
  }
  matrice/nb_probs
}  

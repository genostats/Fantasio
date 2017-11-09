HBD.recap <- function(h, by_segments=F)
{
  proba <- submap.HBD(h)#recuperer les probas HBD pour chaque sous cartes sous forme d'une liste

  if(by_segments)
  {
    nom <- unique(unlist(sapply(proba, function(x) rownames(x)))) #recuperer le nom de chaque individus apparut
    matrice <- matrix(NA, nrow=length(nom), ncol=ncol(proba[[1]]))#creer une matrice avec les individus en lignes 
                                                                  #et les segments en colonnes
    rownames(matrice) <- nom
    Sum <- 0
    for(j in 1 : length(nom))
    {
      cpt <- 0
      for(i in 1:length(proba))
      {
        line <- which(nom[j] == rownames(proba[[i]]))
        if(length(line) == 0) next()
        Sum <- Sum + proba[[i]][line, ]
        cpt <- cpt + 1 
      }
      matrice[j, ] <- Sum / cpt
      
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
  
  po <- match(colnames(matrice), x@snps$id)
  snp.chr <- x@snps$chr[po]
  snp.pos <- x@snps$pos[po]
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

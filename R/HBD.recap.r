HBD.recap <- function(h)
{
  proba <- submap.HBD(h)

  #tableau comptant le nombre de fois ou un marqueur a ete selectionne
  marqueurs <- c()
  for(i in 1:5)
     marqueurs <- c(marqueurs, colnames(proba[[i]]))
  marqueurs <- as.data.frame(table(marqueurs), stringsAsFactors=FALSE)
  
  #tableau comptant le nombre de fois ou un individus est apparu
  individuals <- c()
  for(i in 1:5)
      individuals <- c(individuals, rownames(proba[[i]]))
  individuals <- as.data.frame(table(individuals), stringsAsFactors=FALSE) 
  
  #creation de la matrice = individus x marqueurs
  matrice <- matrix(0, nrow = nrow(individuals), ncol = nrow(marqueurs))
  dimnames(matrice) <- list(individuals[,1], marqueurs[,1])
  
  #trouver l'/les indice(s) et la/les sous-carte(s) avec des valeurs pour chaque marqueur
  az <- vapply(proba, function(jj) match(marqueurs[,1], colnames(jj)), integer(length(marqueurs[,1])))
  
  #boucle de remplissage de la matrice
  for( i in 1:ncol(matrice))
  {
    az1 <- az[i,]
    I <- which(!is.na(az1))
    v <- sapply(I, function(i) proba[[i]][, az1[i]], simplify = F)
    
    #if(is.list(v))
    #{
    division <- c()# conserver tous les indices des individus
      for( j in 1:length(v))
      {
        m <- match(rownames(as.data.frame(v[[j]])), rownames(matrice))
        division <- c(division, m)
        matrice[m,i] <- matrice[m,i] + v[[j]]
      }
    matrice[division,i] <- matrice[division,i] / marqueurs[i,2] #diviser par le nombre de fois ou le marqueur a ete selectionne
  }
    #else{
    #  m <- match(rownames(v), rownames(matrice))
    #  matrice[m,i] <- v
    #}
    matrice
  }  

read.msat.matrix <- function(mapfile, datafile) {
  
  # premiC(re lecture de mapfile -> compter nb lignes et nb colonnes du fichier mapfile
  
  a <- file(mapfile, open = "r") # ouverture du fichier mapfile 
  
  maxcol <- 0 #compteur nb de colonne du fichier
  maxrow <- 0 #compteur nb de ligne du fichier
  
  while(TRUE)
  {
    L <- scan(a, what = "character", nlines = 1, quiet = TRUE)
    le <- length(L)
    if(le == 0) break;
    maxrow = maxrow + 1
    if(le > maxcol) maxcol <- le
  }
  
  close(a)
  
  # deuxiC(me lecture : crC)ation des matrice des frC)quences et de la matrice des noms des allC(les
  
  freq <- matrix(0.0, nrow = maxrow, ncol = (maxcol - 4)/2) # matrice des frC)quences = nb msat * nb max d'allC(les
  geno <- matrix(" ", nrow = maxrow, ncol = (maxcol - 4)/2) # matrice des noms des allC(les
  map <- data.frame(id = rep(NA, maxrow), chr = 0L, distance = 0.)# crC)ation de la matrice map (contient les 3 premiC(res colonnes)
  
  a <- file(mapfile, open = "r") #fermeture et ouverture pour reprendre le fichier au dC)but
  
  for(i in 1:maxrow)
  {
    L <- scan(a, what = "character", nlines = 1, quiet = TRUE)
    le <- length(L)   
    co <- (le-4) / 2
    freq[i, 1:co] <- as.numeric( L[seq(6, le, by = 2)] ) #remplissage de la matrice des frC)quance
    geno[i, 1:co] <- L[ seq(5, le, by = 2) ] #remplissage de la matrice des noms des allC(les
    map$id[i] <- L[1]
    map$chr[i] <- as.integer(L[2])
    map$distance[i] <- as.numeric(L[3])
  }
  close(a)
  
  #compter le nombre de ligne dans le fichier datafile pour crC)er la matrice des gC)notypes
  b <- file(datafile, open="r") #ouverture du fichier datafile 
  
  maxrow_datafile <- 0 #compter le nb de ligne ds le fichier datafile (voir plus bas pourquoi) *
  
  while(TRUE)
  {
    
    L <- scan(b, what = "character", nlines = 1, quiet = TRUE)
    le <- length(L)
    if(le == 0) break
    maxrow_datafile = maxrow_datafile + 1
    
  }
  
  close(b)
  
  # matrice des gC)notypes = nb de msat(ligne dans le fichier mapfile) * 2(nb d'individus)
  # soit nb de ligne ds le fichier data *
  
  matrice_genotype <- matrix("", nrow = maxrow, ncol = (2 * maxrow_datafile) ) #msat matrix non re numC)rotC)
  ped <- data.frame( famid = rep(NA, maxrow_datafile), 
                     id = NA, father = NA, mother = NA, sex = NA, phenotype = NA) #crC)ation de la matrice ped
  
  b <- file(datafile, open="r")#reprendre le fichier au dC)but
  
  #OBJECTIF
  #lire en ligne le fichier et remplir en colonne la matrice des gC)notypes
  
  for (i in 1:maxrow_datafile)
  {
    
    L <- scan(b, what = "character", nlines = 1, quiet = TRUE)
    le <- length(L)
    if(le == 0) break
    matrice_genotype[, 2*i-1] <- L[seq(7,le, by = 2)]
    matrice_genotype[, 2*i] <- L[seq(8,le, by = 2)]
    ped[i, 1:6] <- L[seq(1,6)]
  }
  
  
  close(b)
  
  #crC)ation de la msat matrix re numC)rotC) grC"ce aux deux matrices crC)er prC)cC)demment 
  #(matrice des gC)notypes/nom des allC(les)
  #les allC(les non prC)sents dans le mapfile sont laissC)s C  0 -> genotype inconnu/manquant
  
  matrice_genotype_modif <- matrix(0L, nrow = maxrow, ncol = (2* maxrow_datafile)) 
  
  for (i in 1:ncol(geno))
  {
    for(j in 1:nrow(geno))
    {
      matrice_genotype_modif [ j, matrice_genotype[j,] == geno[j,i] ] <- i
    }
  }
  
  
  new("msat.matrix", maxrow, maxrow_datafile, ped, matrice_genotype_modif, map, freq)
}



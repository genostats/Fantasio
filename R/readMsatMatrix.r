readMsatMatrix <- function(mapfile, datafile) {
  
  # first read of the mapfile -> count the number of lines and the number of columns of the mapfile.
  
  a <- file(mapfile, open = "r") #opening of the mapfile
  
  maxcol <- 0 #count the number of columns of the file
  maxrow <- 0 #count the number of lines of the file
  
  while(TRUE)
  {
    L <- scan(a, what = "character", nlines = 1, quiet = TRUE)
    le <- length(L)
    if(le == 0) break;
    maxrow = maxrow + 1
    if(le > maxcol) maxcol <- le
  }
  
  close(a)
  
  #second read of the file : creation of frequencies and allelic name matrix
  
  freq <- matrix(0.0, nrow = maxrow, ncol = (maxcol - 4)/2) # frequencies matrix = nb msat * max nb of allele
  geno <- matrix(" ", nrow = maxrow, ncol = (maxcol - 4)/2) # name of allele matrix
  map <- data.frame(id = rep(NA, maxrow), chr = 0L, distance = 0.)# creation of the matrice map, first 3 columns
  
  a <- file(mapfile, open = "r") #closing and opening the file to start from the beginning
  
  for(i in seq_len(maxrow))
  {
    L <- scan(a, what = "character", nlines = 1, quiet = TRUE)
    le <- length(L)   
    co <- (le-4) / 2
    freq[i, seq_len(co)] <- as.numeric( L[seq(6, le, by = 2)] ) #filling the frequencies matrix
    geno[i, seq_len(co)] <- L[ seq(5, le, by = 2) ] #filling the name of allele matrix
    map$id[i] <- L[1]
    map$chr[i] <- as.integer(L[2])
    map$distance[i] <- as.numeric(L[3])
  }
  close(a)
  
  #count the number of lines on the datafile to create the genotype matrix
  b <- file(datafile, open="r") #open datafile
  
  maxrow_datafile <- 0 #count the number of lines on datafile
  
  while(TRUE)
  {
    
    L <- scan(b, what = "character", nlines = 1, quiet = TRUE)
    le <- length(L)
    if(le == 0) break
    maxrow_datafile = maxrow_datafile + 1
    
  }
  
  close(b)
  
  # genotype matrix = nb msat * 2 nb of individual
  
  matrice_genotype <- matrix("", nrow = maxrow, ncol = (2 * maxrow_datafile) ) #original msat matrix 
  ped <- data.frame( famid = rep(NA, maxrow_datafile), 
                     id = NA, father = NA, mother = NA, sex = NA, phenotype = NA) #ped matrix
  
  b <- file(datafile, open="r")
  
  #OBJECTIV
  #read by line the file et fill by columns the genotype matrix
  
  for (i in seq_len(maxrow_datafile))
  {
    
    L <- scan(b, what = "character", nlines = 1, quiet = TRUE)
    le <- length(L)
    if(le == 0) break
    matrice_genotype[, 2*i-1] <- L[seq(7,le, by = 2)]
    matrice_genotype[, 2*i] <- L[seq(8,le, by = 2)]
    ped[i, seq_len(6)] <- L[seq(1,6)]
  }
  
  
  close(b)
  
  #creation of the msat matrix numeroted thanks to previous created matrix 
  #non present allel are left to 0 <- unknown/missing genotype
  
  matrice_genotype_modif <- matrix(0L, nrow = maxrow, ncol = (2* maxrow_datafile)) 
  
  for (i in seq_len(ncol(geno)))
  {
    for(j in seq_len(nrow(geno)))
    {
      matrice_genotype_modif [ j, matrice_genotype[j,] == geno[j,i] ] <- i
    }
  }
  
  
  new("msat.matrix", maxrow, maxrow_datafile, ped, matrice_genotype_modif, map, freq)
}



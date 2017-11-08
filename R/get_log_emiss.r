# x est une msat_matrix
get.log.emiss <- function(x, i) #lignes de log.emiss qui concernent l'individu i
{
  x@log.emiss[ c(2*i-1, 2*i), ]
}
  

##################################################################################
#This function create a summary for HFLOD score                                  #
#                                                                                #
#!!! submaps : the list object                                                   #                                       
#                                                                                #
#*** return a matrix with HFLOD score                                            #
##################################################################################

#
#HFLOD.recap <- function(submaps)
#{
#  #if one submap only return it
#  if(length(submaps@atlas) == 1)
#  {
#    matrice <- submaps@HFLOD
#    return(matrice)	
#  }
#  
#  #tous les scores HFLOD sur toutes les sous-cartes
#  
#  matrice <- matrix(NA, nrow = nrow(HFLOD_scores[[1]]),
#                        ncol = ncol(HFLOD_scores[[1]]))
#  
#  
#  #moyenne des HFLOD scores et de moving average par segment
#  for( i in 1:nrow(matrice))
#  {
#    HFLOD      <- mean(sapply(HFLOD_scores, function(hh) hh[i,1]))
#    moving_Avg <- mean(sapply(HFLOD_scores, function(hh) hh[i,2]))
#    matrice[i, 1] <- HFLOD
#    matrice[i, 2] <- moving_Avg
#  }
#  
#  return(matrice)
#}#
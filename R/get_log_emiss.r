
##################################################################################
#This function rearrange the log emission to get a dataframe                     #
# 2 x the_of_individuals-1 & 2xthe_of_individuals                                #
#                                                                                #
#!!! x : a submap matrix                                                         #                                       
#!!! i : the individual i                                                        #
#                                                                                #
#*** return the submap matrix with the modified dataframe                        #
##################################################################################

get.log.emiss <- function(x, i) #lignes de log.emiss qui concernent l'individu i
{
  x@log.emiss[ c(2*i-1, 2*i), ]
}
  

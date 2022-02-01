
##################################################################################
#This function rearrange the log emission to get a dataframe                     #
# 2 x the_of_individuals-1 & 2xthe_of_individuals                                #
#                                                                                #
#!!! x : a submap matrix                                                         #                                       
#!!! i : the individual i                                                        #
#                                                                                #
#*** return the submap matrix with the modified dataframe                        #
##################################################################################

getLogEmiss <- function(x, i) #line of log.emmis for the individual i 
{
  x@log.emiss[ c(2*i-1, 2*i), ]
}
  

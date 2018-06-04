
##################################################################################
#This function computes the log.emission for a submaps                           #
#by calling a function written in C++                                           #
#                                                                                #
#!!! x : a bedmatrix                                                             #                                       
#!!! map : a submap matrix object                                                #
#!!! epsilon : the value of epsilon for submap                                   #
#                                                                                #
#*** a dataframe with log.emiss computed                                         #
##################################################################################

bed.logEmiss <- function(x, map, epsilon = 1e-3) {
  
  if(epsilon < 0 | epsilon > 1)
    stop("The genotyping error rate 'epsilon' should be between 0 and 1")
  if(epsilon > 0.01)
    warning("A high genotyping error rate 'epsilon' is unlikely to give sensible results")
  .Call('festim_m4_logEmiss', PACKAGE = "FEstim", x@bed, x@p, map, epsilon)
  
}
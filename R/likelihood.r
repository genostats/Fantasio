##################################################################################
#This function computed likelihood0 and 1 for a marker                           #
#                                                                                #
#!!! logEmission : the value of logEmission for the marker                       #                                       
#!!! delta.Distance : the value of delta.distance for the marker                 #
#!!! a : the value of a for the marker                                           #
#!!! f : the value of f for the marker                                           #
#                                                                                #
#*** return the value of likelihood                                            #
##################################################################################

Likelihood <- function(logEmission, delta.Distance, a, f) {
  if( a < 0 ) stop("a should be >= 0")
  if( f < 0 | f > 1) stop("f should be between 0 and 1")
  .Call('festim_logLikelihood', PACKAGE = "Fantasio", logEmission, delta.Distance, a, f)
}

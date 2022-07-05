##################################################################################
#This function computes HBD probabiilties for a marker                           #
#                                                                                #
#!!! logEmission : the logEmission for the marker                                #                                       
#!!! delta.Distance : the value of delta.distance for the maker                  #
#!!! a : the value of a for the marker                                           #
#!!! f : the value of f for the marker                                           #
#                                                                                #
#*** return the value of HBD                                                     #
##################################################################################

likelihoodGradient <- function(logEmission, delta.Distance, a, f) {
  if( a < 0 ) stop("a should be >= 0")
  if( f < 0 | f > 1) stop("f should be between 0 and 1")
  .Call('_Fantasio_logLikelihood_gradient', PACKAGE = "Fantasio", logEmission, delta.Distance, a, f)
}


log.proba.emission <- function(Y1, Y2, logFreq, epsilon = 1e-3) {
  N <- length(Y1)
  if(length(Y2) != N) stop()
  if(nrow(logFreq) != N) stop()
  if(max(Y1, Y2) > ncol(logFreq)) stop()

  .Call('festim_logEmiss', PACKAGE = "Fantasio", Y1, Y2, logFreq, epsilon)
}

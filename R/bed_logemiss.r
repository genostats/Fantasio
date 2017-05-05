
bed.logEmiss <- function(x, map, epsilon = 1e-3) {
  .Call('festim_m4_logEmiss', PACKAGE = "FEstim", x@bed, x@p, map, epsilon)
}

#' @import gaston
#' @import methods
#' @import parallel
#' 
#' @importFrom graphics abline axis lines par plot points polygon text rect
#' @importFrom grDevices dev.off png gray
#' @importFrom methods callNextMethod new
#' @importFrom stats median optim optimize pchisq 
#' @importFrom utils data read.table tail
#' 
#' @useDynLib Fantasio
#' 
#' @exportClass f.matrix
#' @exportClass list.submaps
#' @exportClass hotspots.matrix
#' @exportClass hotspot.segments
#' @exportClass snps.matrix
#' @exportClass snps.segments

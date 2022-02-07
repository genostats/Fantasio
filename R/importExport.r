#' @useDynLib Fantasio
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
#' 
#' @exportClass fMatrix
#' @exportClass atlas
#' @exportClass HostspotsMatrix
#' @exportClass HostspotsSegments
#' @exportClass snpsMatrix
#' @exportClass snpsSegments
NULL
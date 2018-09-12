setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
#' Class fMatrix.
#'
#' Class defines a mother class for all the subclasses in the package Fantasio.
#'
#' @rdname fMatrix-class
#' @exportClass fMatrix
#' @slot ncol number of loci
#' @slot nrow number of individual
#' @slot ped  the first 6 columns of a .ped
#' @slot map  id, chr, distance
#' @slot epsilon value of epsilon = genotyping error rare, to use for computing log emission precalculated at initialisation
#' @slot delta.dist difference of map$distance 
#' @slot log.emiss matrix of log proba emission if statut = 0 or 1(2 nb inds x nb msats) 
#' @slot a value of a estimated by festim
#' @slot f value of f estimated by festim
#' @slot likelihood0 likelihood under H0 (f = 0)
#' @slot likelihood1 likelihood under H1 
#' @slot p.lrt likelihood ratio test
#' @slot HBD.prob proba HBD = 1 ; one individual per column : dim = (nb inds x nb msats)
#' @slot FLOD matrix of FLOD scores dim = (nb inds x nb msats)
setClass("fMatrix", representation(
          # input slots
          ncol = 'numeric',              # number of loci
          nrow = 'numeric',              # number of individual
          ped = 'data.frame',            # the first 6 columns of a .ped
          map = 'data.frame',            # id, chr, distance
          epsilon = 'numeric',           # value of epsilon = genotyping error rare, to use for computing log emission precalculated at initialisation
          delta.dist = 'numericOrNULL',  # diff(map$distance) + faire attention au chgt de chr
          log.emiss = 'matrixOrNULL',    # matrix of log proba emission if statut = 0 or 1(2 nb inds x nb msats) 
          a = 'numeric',                 # value of a and f estimated by festim
          f = 'numeric',
          likelihood0 = 'numeric',       # likelihood under H0 (f = 0)
          likelihood1 = 'numeric',       # likelihood under H1 
          p.lrt = 'numeric',             # likelihood ratio test
          HBD.prob = 'matrix',           # proba HBD = 1 ; one individual per column : dim = (nb inds x nb msats)
          FLOD = 'matrix'                # matrix of FLOD scores dim = (nb inds x nb msats)
))

#' Constructor method of fMatrix
#'
#' @param .Object the object type
#' @param ncol number of loci
#' @param nrow number of individual
#' @param ped  the first 6 columns of a .ped
#' @param map  id, chr, distance
#' @param log.emiss matrix of log proba emission if statut = 0 or 1(2 nb inds x nb msats) 
#' @param epsilon value of epsilon = genotyping error rare, to use for computing log emission precalculated at initialisation
setMethod('initialize', signature='fMatrix', definition=function(.Object, ncol, nrow, ped, map, log.emiss, epsilon=1e-3) {
  if(nrow(ped) != nrow) stop("ped dimension mismatch")
  if(nrow(map) != ncol) stop("map dimension mismatch")

  delta.dist <- diff(map$distance)
  # properly treat the change of chromosome
  I <- cumsum(rle(map$chr)$length)
  I <- I[-length(I)]
  delta.dist[I] <- -1 
  
  .Object@ncol     <- ncol
  .Object@nrow     <- nrow
  .Object@ped      <- ped
  .Object@map      <- map 
  .Object@delta.dist <- delta.dist
  .Object@log.emiss <- log.emiss
  .Object@epsilon   <- epsilon
  .Object
})

setGeneric('dim')

#' The generic dim method for fMatrix class
#'
#' @param x a bed.matrix object
setMethod("dim", signature = "fMatrix", function(x) c(x@nrow, x@ncol) )

#' show method of fMatrix class
#' 
#' @param object an fMatrix object
setMethod('show', signature("fMatrix"),
          function(object) {
            cat('A fMatrix with ', nrow(object), ' individual(s) and ', ncol(object), " markers\n")
          }
)

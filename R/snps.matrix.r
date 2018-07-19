#' Class snps.matrix
#'
#' Class \code{snps.matrix} This class is use to create an object that will represent the submap created    
#'
#' @rdname snps.matrix-class
#' @exportClass snps.matrix
#' @slot step the value of the step used to created the segments
#' @slot submap a vector of index
setClass("snps.matrix", representation(
  step = "numeric",
  submap = 'numeric'              # vecteur d'indices
), contains = "f.matrix" )

#' Constructor method of snps.matrix.
#'
#' @param .Object the object type
#' @param step the step use to take step in a mini segments (refer to vignette for more infos)
#' @param ncol number of loci
#' @param nrow number of individual
#' @param submap a list of submaps
#' @param ped  the first 6 columns of a .ped
#' @param map  id, chr, distance
#' @param log.emiss matrix of log proba emission if statut = 0 or 1(2 nb inds x nb msats) 
#' @param epsilon value of epsilon = genotyping error rare, to use for computing log emission precalculated at initialisation
setMethod('initialize', signature='snps.matrix', definition=function(.Object, step, ncol, nrow, submap, ped, map, log.emiss, epsilon = 1e-3) {
  .Object@step <- step
  .Object@submap  <-  submap
  callNextMethod(.Object, ncol, nrow, ped, map, log.emiss, epsilon)
})

#' Show method of snps.matrix.
#'
#' @param object an snps.matrix object
setMethod('show', signature("snps.matrix"),
          function(object) {
            cat('A snps.matrix with ', nrow(object), ' individual(s) and ', ncol(object), " markers\n")
          }
)



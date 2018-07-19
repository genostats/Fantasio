#' Class snps.segments.
#'
#' Class \code{snps.segments} defines snps.segments object that will store the list of segments created.
#'
#' @rdname snps.segments-class
#' @exportClass snps.segments
#' @slot gap the value of the gap used to created the segments
#' @slot unit the unit of the marker in the segments (cM or bP)
#' @slot snps.segments the list of segments
setClass("snps.segments", representation(
  gap = 'numeric',
  unit = 'character',
  snps.segments = 'list'
))

#' Constructor method of snps.segments.
#'
#' @param .Object the type of the object ( here snps.segments)
#' @param gap the gap used to creates the segments 
#' @param unit the unit in which are the markers
#' @param segments the list of segments
setMethod('initialize', signature='snps.segments', definition=function(.Object, gap, unit, segments) {
  .Object@gap <- gap
  .Object@unit <- unit
  .Object@snps.segments <- segments
  .Object
})

#' Show method of snps.segments.
#'
#' @param object the snps.segments object
setMethod('show', signature("snps.segments"), 
          function(object){
            cat('A  list of ', length(object@snps.segments), 'chromosoms segments created using gap in the genome\n')
          })























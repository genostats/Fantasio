#' Class snpsSegments.
#'
#' Class \code{snpsSegments} defines snpsSegments object that will store the list of segments created.
#'
#' @rdname snpsSegments-class
#' @exportClass snpsSegments
#' @slot gap the value of the gap used to created the segments
#' @slot unit the unit of the marker in the segments (cM or bP)
#' @slot snpsSegments the list of segments
setClass("snpsSegments", representation(
  gap = 'numeric',
  unit = 'character',
  snpsSegments = 'list'
))

#' Constructor method of snpsSegments.
#'
#' @param .Object the type of the object ( here snpsSegments)
#' @param gap the gap used to creates the segments 
#' @param unit the unit in which are the markers
#' @param segments the list of segments
setMethod('initialize', signature='snpsSegments', definition=function(.Object, gap, unit, segments) {
  .Object@gap <- gap
  .Object@unit <- unit
  .Object@snpsSegments <- segments
  .Object
})

#' Show method of snpsSegments.
#'
#' @param object the snpsSegments object
setMethod('show', signature("snpsSegments"), 
          function(object){
            cat('A  list of ', length(object@snpsSegments), 'chromosoms segments created using gap in the genome\n')
          })























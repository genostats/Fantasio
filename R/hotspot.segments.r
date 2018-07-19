##################################################################################
#This is the class that encapsulate the list of segments created by using the    #
#hotspots in the genome                                                          #
#                                                                                #
##################################################################################
#' Class hotspot.segments
#'
#' Class \code{hotspot.segments} This is the class that encapsulate the list of segments created by using the hotspots in the genome.
#' @exportClass hotspot.segments
#' 
setClass("hotspot.segments", representation = (
  hotspot.segments = 'list'
))

#' Show method of hotspot.segments.
#'
#' @param object the hotspot.segments
setMethod('show', signature("hotspot.segments"), 
          function(object){
            cat('A  list of ', length(object), 'chromosoms segments created using hotspots in the genome\n')
          })























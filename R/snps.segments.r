##################################################################################
#This is the class that encapsulate the list of segments created by using the gap#
#in the genome                                                                   #
#                                                                                #
##################################################################################

setClass("snps.segments", representation(
  gap = 'numeric',
  snps.segments = 'list'
))

setMethod('initialize', signature='snps.segments', definition=function(.Object, gap, segments) {
  .Object@gap <- gap
  .Object@snps.segments <- segments
  .Object
})

setMethod('show', signature("snps.segments"), 
          function(object){
            cat('A  list of ', length(object@snps.segments), 'chromosoms segments created using gap in the genome\n')
          })























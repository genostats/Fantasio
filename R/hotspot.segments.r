##################################################################################
#This is the class that encapsulate the list of segments created by using the    #
#hotspots in the genome                                                          #
#                                                                                #
##################################################################################

setClass("hotspot.segments", representation = (
  hotspot.segments = 'list'
))
setMethod('show', signature("hotspot.segments"), 
          function(object){
            cat('A  list of ', length(object), 'chromosoms segments created using hotspots in the genome\n')
          })























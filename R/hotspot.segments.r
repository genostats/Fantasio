setClass("hotspot.segments", representation = (
  hotspot.segments = 'list'
))
setMethod('show', signature("hotspot.segments"), 
          function(object){
            cat('A  list of ', length(object), 'chromosoms, each slot is a list which contains index of each markers between two hotspots (ie: in a segment)\n')
          })






























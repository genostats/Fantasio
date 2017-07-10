setClass("list.submaps", representation = (
        submap.info = 'list'
))
setMethod('show', signature("list.submaps"), 
  function(object){
       cat('A list of submaps with ', length(object), ' submaps\n ')
  })




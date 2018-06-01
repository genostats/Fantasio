setClass("submap.matrix", representation(
       submap = 'numeric'              # vecteur d'indices
), contains = "f.matrix" )

setMethod('initialize', signature='submap.matrix', definition=function(.Object, ncol, nrow, submap, ped, map, log.emiss, epsilon = 1e-3) {
  .Object@submap  <-  submap
  callNextMethod(.Object, ncol, nrow, ped, map, log.emiss, epsilon)
})

setMethod('show', signature("submap.matrix"),
  function(object) {
    cat('A submap.matrix with ', nrow(object), ' individual(s) and ', ncol(object), " markers\n")
  }
)



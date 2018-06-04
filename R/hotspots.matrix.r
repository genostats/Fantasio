##################################################################################
#This function is the class that is use to create an object that will represent  #
#the submap created                                                              #
#                                                                                #
##################################################################################

setClass("hotspots.matrix", representation(
       submap = 'numeric'              # vecteur d'indices
), contains = "f.matrix" )

setMethod('initialize', signature='hotspots.matrix', definition=function(.Object, ncol, nrow, submap, ped, map, log.emiss, epsilon = 1e-3) {
  .Object@submap  <-  submap
  callNextMethod(.Object, ncol, nrow, ped, map, log.emiss, epsilon)
})

setMethod('show', signature("hotspots.matrix"), function(object) {
    cat('A hotspots.matrix with ', nrow(object), ' individual(s) and ', ncol(object), " markers\n")
  }
)



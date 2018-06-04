##################################################################################
#This function is the class that is use to create an object that will represent  #
#the submap created                                                              #
#                                                                                #
##################################################################################

setClass("snps.matrix", representation(
  step = "numeric",
  submap = 'numeric'              # vecteur d'indices
), contains = "f.matrix" )

setMethod('initialize', signature='snps.matrix', definition=function(.Object, step, ncol, nrow, submap, ped, map, log.emiss, epsilon = 1e-3) {
  .Object@step <- step
  .Object@submap  <-  submap
  callNextMethod(.Object, ncol, nrow, ped, map, log.emiss, epsilon)
})

setMethod('show', signature("snps.matrix"),
          function(object) {
            cat('A snps.matrix with ', nrow(object), ' individual(s) and ', ncol(object), " markers\n")
          }
)



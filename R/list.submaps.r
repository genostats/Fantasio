#######################################################################################################################
#This is the class used to create an object which will contains every dataframe and list created when creating submaps#
#######################################################################################################################

setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("dataframeOrNULL",members=c("data.frame", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("characterOrNULL",members = c("character", "NULL"))
setClassUnion("doubleOrNULL",members = c("numeric", "NULL"))


setClass("list.submaps", representation(
        segments_list        = 'listOrNULL',
        atlas                = 'list', 
        likelihood_summary   = 'listOrNULL',
        estimation_summary   = 'listOrNULL', 
        marker_summary       = 'dataframeOrNULL',
        submap_summary       = 'dataframeOrNULL',
        HBD_recap            = 'matrixOrNULL',
        FLOD_recap           = 'matrixOrNULL',  
        HBD_segments         = 'listOrNULL',
        HFLOD                = 'dataframeOrNULL',
        bedmatrix            = 'bed.matrix', 
        bySegments           = "logical", 
        unit                 = "characterOrNULL", 
        gap                  = "doubleOrNULL"
))
setMethod('initialize', signature='list.submaps', definition=function(.Object, submaps, bedmatrix, segments_list, bySegments, unit=NULL, gap=NULL)
{
  .Object@atlas         <- submaps
  .Object@bedmatrix     <- bedmatrix
  .Object@segments_list <- segments_list
  .Object@bySegments    <- bySegments
  .Object@unit          <- unit
  .Object@gap           <- gap
  .Object
})
setMethod('show', signature("list.submaps"), 
  function(object){
       cat('A list of', length(object@atlas), 'submaps\n ')
  })




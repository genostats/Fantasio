setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("dataframeOrNULL",members=c("data.frame", "NULL"))
setClass("list.submaps", representation(
        atlas = 'list', 
        likelihood_summary = 'listOrNULL',
        estimation_summary = 'listOrNULL', 
        marker_summary = 'dataframeOrNULL',
        submap_summary = 'dataframeOrNULL'
))
setMethod('initialize', signature='list.submaps', definition=function(.Object, submaps)
{
  .Object@atlas <- submaps
  .Object
})
setMethod('show', signature("list.submaps"), 
  function(object){
       cat('A list of', length(object@atlas), 'submaps\n ')
  })




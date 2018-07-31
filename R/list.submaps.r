#######################################################################################################################
#This is the class used to create an object which will contains every dataframe and list created when creating submaps#
#######################################################################################################################

setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("dataframeOrNULL",members=c("data.frame", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("characterOrNULL",members = c("character", "NULL"))
setClassUnion("doubleOrNULL",members = c("numeric", "NULL"))

#' Class list.submaps
#'
#' Class \code{list.submaps} This is the class used to create an object which will contains every dataframe and list created when creating submaps.  
#'
#' @rdname list.submaps-class
#' @exportClass list.submaps
#' @slot segments_list the list of segments 
#' @slot atlas a list of submaps
#' @slot likelihood_summary a dataframe with both likelihood0 and likelihood1 over the submaps
#' @slot estimation_summary a dataframe with both a and f estimation over the submaps
#' @slot marker_summary a dataframe indicating the number of times a marker has been chosen
#' @slot submap_summary a dataframe with summary statitistics about the submaps
#' @slot HBD_recap a dataframe with for one individual and for one marker a mean computation of all the HBD probabilities computed, on every individuals.
#' @slot FLOD_recap a dataframe with for one individual and for one marker a mean computation of all the FLOD scores computed, on every individuals.
#' @slot HBD_segments a list of dataframe with the HBD segments, on all individuals.
#' @slot HFLOD a dataframe with the value of HFLOD scores for every markers through all submaps.
#' @slot bedmatrix  a bed.matrix object (refer to gaston package)
#' @slot bySegments a boolean indicating wheter the creation of summary statistics was made by segments (see documentation of Fantasio function)
#' @slot unit   the unit of the markers (cM or Bp).
#' @slot gap   the value of the gap used to pick marker when doing submaps by snps. (see function Fantasio for more infos)
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

#' Constructor method of list.submaps.
#'
#' @name list.submaps-class
#' @rdname list.submaps-class
#' @param .Object the object type
#' @param submaps a list of submaps
#' @param bedmatrix a bed.matrix object
#' @param segments_list a list of segments object
#' @param bySegments a boolean
#' @param unit the unit in which are markers (cM or Bp)
#' @param gap the gap used to create segments in the By Distance method
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


#' Show method of list.submaps.
#'
#' @param object an list.submaps object
setMethod('show', signature("list.submaps"), 
  function(object){
       cat('A list of', length(object@atlas), 'submaps\n ')
  })




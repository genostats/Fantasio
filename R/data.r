#' Fantasio.
#'
#' @name Fantasio
#' @docType package
#' 
NULL

#' Position of hotspots throughout the human genome 
#'
#' These datasets are data frames containing position and intensity of hotspots in the human genome,
#' for releases hg17 (build 35), hg18 (build 36) and hg19 (build 37).
#'
#' The variables are:
#' \itemize{
#'   \item Chromosome. Chromosome number
#'   \item Centre. Center of the hotspots in Mb
#'   \item Start. Start of the hotspot
#'   \item End.  End of the hotspot
#'   \item Widthkb. Width of the hotposts in kb
#'   \item IntensitycMMb. Intensity in M / Mb
#'   \item TotaldistancecM. Total distance in cM
#' }
#'
#' @docType data
#' @keywords datasets
#' @name hotspot
#' @aliases hotspot_hg17 hotspot_hg18 hotspot_hg19
#' @usage data(hotspot_hg17)
#' @format Data frames with 32996 rows and 7 variables
NULL

#' Position of cytobands on the human chromosomes
#'
#' These datasets are data frames containing position of cytobands on human
#' chromosomes, for builds 35 (hg17), 36 (hg18), 37 (hg19) and 38 respectively.
#' The data were build starting from the data given in Bioconductor package 
#' \code{quantsmooth}.
#'
#' The variables are:
#' \itemize{
#'   \item chr     Chromosome number.
#'   \item arm     Chromosome arm.
#'   \item band    Band name.
#'   \item stain   Band stain.
#'   \item beg.bp  Beginning position in base pairs.
#'   \item end.bp  Ending position in base pairs.
#'   \item beg.cM  Beginning position in centiMorgans.
#'   \item end.cM  Ending position in centiMorgans.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name cytobands
#' @aliases cytobands.b35  cytobands.b36  cytobands.b37  cytobands.b38  
#' @usage data(cytobands.b35)
#' @format A data frame with 32996 rows and 7 variables
NULL

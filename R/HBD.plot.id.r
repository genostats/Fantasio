#' plot of HBD segment 
#' 
#' This function plots the HBD segments for all the chromosoms of a given individual
#' 
#' @param Submaps a list.submap object
#' @param ROHfile a ROH file from which the segments will be plotted (optional)
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
#' @param individual.id the individual id of the individual wanted
#' @param family.id the family id of the individual wanted
#' @param regions a specific region to be enlighted in the plot (optional)
#' @param outfile a name for the plot (optional)
#' @param build the value of the build to use to plot chromosome in the plot value accepted are 35, 36, 37, 38 (default is 37)
#' 
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end
#' 
#' @return return a plot of the individual's HBD segments.
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' 
#' @export
HBD.plot.id <- function(Submaps, ROHfile, unit= "cM", individual.id, family.id, regions, outfile, build = 37)
{
  if(class(Submaps@atlas[[1]])[1] != "snps.matrix" & class(Submaps@atlas[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps.") 
  
  if(class(Submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix.")
  
  if(is.null(Submaps@HBD_recap))
    stop("HBD_recap is empty cannot plot, make sure to have atleast one individual considered INBRED.")
  
  if(!missing(Submaps) & !missing(ROHfile))
  {
    plot.ROH.segments.id(Submaps=Submaps, ROHfile, unit, individual.id=individual.id, family.id=family.id, regions, outfile=outfile, build=build)
  }else{
    if(!missing(Submaps))
      plot.HBD.segments.id(Submaps = Submaps, individual.id=individual.id, family.id=family.id, unit=unit, regions, outfile=outfile, build=build)
  }
}
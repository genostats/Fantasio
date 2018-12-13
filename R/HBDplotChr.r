#' plot of HBD segment 
#' 
#' This function plots the HBDsegments for a given chromosom and all the individual
#' 
#' @param Submaps a list.submap object
#' @param ROHfile a ROH file from which the segments will be plotted (optional)
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
#' @param chr the chromosome number from which to plot HBD segment
#' @param list.ids a vector containing a list of individuals from which only the HBDsegments for this chromosome will be ploted (optional)
#' @param regions a specific region to be enlighted in the plot (optional)
#' @param outfile a name for the plot (optional)
#' @param inbred whether you want to plot only INBRED individuals or not (default is FALSE)
#' @param build the value of the build to use to plot chromosome in the plot value accepted are 35, 36, 37, 38 (default is 37)
#' 
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end
#' 
#' @seealso Fantasio
#' @seealso setHFLOD
#' 
#' @return return a plot of the chromosome HBDsegments for all the individual
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' 
#' @export
HBDplotChr <- function(Submaps, ROHfile, unit="cM", chr, list.ids, regions, outfile, inbred = FALSE, build=37)
{
  if(inbred)
  {
    list.ids <- which(Submaps@submap_summary$INBRED ==TRUE)
    ind      <- as.vector(Submaps@submap_summary$IID[list.ids])
    fam      <- as.vector(Submaps@submap_summary$FID[list.ids])
    list.ids <- paste(ind, fam, sep="_")
  }
  
  if(class(Submaps@submaps_list[[1]])[1] != "snpsMatrix" & class(Submaps@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  if(class(Submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix")
  
  if(is.null(Submaps@HBD_recap))
    stop("HBD_recap is empty cannot plot, make sure to have atleast one individual considered INBRED.")
  
  if(!missing(Submaps) & !missing(ROHfile))
  {
    plot.ROH.segments.chr(ROHfile = ROHfile, submaps = Submaps, unit = unit, chr = chr, outfile=outfile, listid=list.ids, regions=regions, build=build)
  }else{
    if(!missing(Submaps))
      plot.HBDsegments.chr(Submaps=Submaps, unit=unit, chr=chr, list.ids=list.ids, regions=regions, outfile=outfile, build=build)
  }
}

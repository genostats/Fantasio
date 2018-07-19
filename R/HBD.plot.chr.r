#' plot of HBD segment 
#' 
#' This function plots the HBD segments for a given chromosom and all the individual
#' 
#' @param Submaps a list.submap object
#' @param ROHfile a ROH file from which the segments will be plotted (optional)
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
#' @param chr the chromosome number from which to plot HBD segment
#' @param list.ids a vector containing a list of individuals from which only the HBD segments for this chromosome will be ploted (optional)
#' @param regions a specific region to be enlighted in the plot (optional)
#' @param outfile a name for the plot (optional)
#' 
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end
#' 
#' @seealso Fantasio
#' @seealso set.HFLOD
#' 
#' @return return a plot of the chromosome HBD segments for all the individual
#' 
#' @examples  
#' ##install.packages("HGDP.CEPH", repos="https://genostats.github.io/R/") ## make this only one time
#' require(Fantasio)
#' require(HGDP.CEPH)
#' filepath <-system.file("extdata", "hgdp_ceph.bed", package="HGDP.CEPH")
#' x <- read.bed.matrix(filepath)
#' x <- set.stats(x)
#' x.me <- select.inds(x, population == "Bedouin")
#' x.me@ped$pheno <- rep(2,48) #The package analyzes only individualw with a status of 2
#' submaps <- Fantasio(x.me, "Hotspots", n=5)
#' HBD.plot.chr(submaps, chr=1)
#' 
#' @export
HBD.plot.chr <- function(Submaps, ROHfile, unit="cM", chr, list.ids, regions, outfile)
{
  
  if(class(Submaps@atlas[[1]])[1] != "snps.matrix" & class(Submaps@atlas[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.") 
  
  if(class(Submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix to eat")
  
  if(is.null(Submaps@HBD_recap))
    stop("HBD_recap is empty cannot plot, make sure to have atleast one individual considered INBRED.")
  
  if(!missing(Submaps) & !missing(ROHfile))
  {
    plot.ROH.segments.chr(ROHfile = ROHfile, submaps = Submaps, unit = unit, chr = chr, outfile=outfile, listid=list.ids, regions=regions)
  }else{
    if(!missing(Submaps))
      plot.HBD.segments.chr(Submaps=Submaps, unit=unit, chr=chr, list.ids=list.ids, regions=regions, outfile=outfile)
  }
}
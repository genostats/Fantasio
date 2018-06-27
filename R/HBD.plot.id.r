#' plot of HBD segment 
#' 
#' This function plots the HBD segments for all the chromosoms of a given individual
#' 
#' @param Submaps a list.submap object
#' @param ROHfile a ROH file from which the segments will be plotted (optional)
#' @param unit the unit uses for position in the plot, "cM" or "Bases"
#' @param individual.id the individual id 
#' @param family.id the family id 
#' @param regions a specific region to be enlighted in the plot (optional)
#' @param outfile a name for the plot (optional)
#' 
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end
#' 
#' @return return a plot of the individual's HBD segments.
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListBySnps(bedMatrix)
#' submaps <- makeSubmapsBySnps(bedMatrix, 5, segmentList)
#' HBD.plot.id(submaps, individual.id="IID", family.id="FID")
#' 
#' @export
HBD.plot.id <- function(Submaps, ROHfile, unit= "cM", individual.id, family.id, regions, outfile)
{
  if(!missing(Submaps) & !missing(ROHfile))
  {
    plot.ROH.segments.id(Submaps=Submaps, ROHfile, unit, individual.id=individual.id, family.id=family.id, regions, outfile=outfile)
  }else{
    if(!missing(Submaps))
      plot.HBD.segments.id(Submaps = Submaps, individual.id=individual.id, family.id=family.id, unit=unit, regions, outfile=outfile)
  }
}

  
  #if(!is.character(individual.id))
  #  return("Need individual.id as character")
  #
  #HBD.recap <- Submaps@HBD_recap
  #HBD.segments <- Submaps@HBD_segments
  #
  ##recuperer les id de l'individus necessaire pour la fonction plot
  #name <- rownames(HBD.recap)
  #id   <- which(name == individual.id)
  #
  ##traitement de l'options regions
  #if (missing(regions)) 
  #  myreg <- NULL
  #else
  #  myreg <- regions
  #
  ##donner un nom au fichier creer
  #if (missing(outfile)) 
  #  outfile <- paste("HBD_", name[id],"_",unit,".png",sep="")
  #else {
  #  outfile <- paste(outfile,".png",sep="") 
  #}
#
  ##creation d'un fichier png et plot
  ##if(save_img)
  ##{
  ##  png(filename = outfile, width = 1000, height = 1000,pointsize=24)
  ##  plot_HBD_IID(HBD.segments, id, unit = unit, regions = myreg, main=paste("HBD segments of ",name[id]))
  ##  dev.off()
  ##}
  #
  #plot_HBD_IID(HBD.segments, id, unit = unit, regions = myreg, main=paste("HBD segments of ",name[id]))

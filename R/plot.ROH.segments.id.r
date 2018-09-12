#' Creation of plot using ROH file
#' 
#' This function is use to plot HBDsegments for a specific infdividual using an ROH file
#' 
#' @param ROHfile an ROH file
#' @param Submaps a submapsList object
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
#' @param regions a matrix containing the value to ve highlighted in the plot
#' @param outfile the name of the plot
#' @param family.id the individual id of the individual wanted
#' @param individual.id the family id of the individual wanted
#' @param build the value of the build to use to plot chromosome in the plot value accepted are 35, 36, 37, 38 (default is 37)
#' 
#' @details Use this function when you want to use an ROH file to plot your HBDsegments for a specific individual.
#' @details Two unit are accepted : "Bases" or "cM".
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end
#'  
#' @return This function returns a plot of the HBDsegments for a specific individual
#' 
plot.ROH.segments.id <- function(Submaps, ROHfile, unit="cM", regions, outfile, family.id, individual.id, build)
{
  ROH <- read.table(ROHfile,header=TRUE)
  ROH <- subset(ROH,ROH$FID==family.id & ROH$IID==individual.id)
  
  if(nrow(ROH) == 0)
    return("No information found for this particular individual, please enter a right family.id and individual.id, and make sure individua_id is a character string")
  
  if (missing(regions))
    myreg <- NULL
  else{
    #myreg <- regions[regions$chr == chr,]
    myreg <- regions
  }
  
  if (unit=="cM")
    ROH <- add_cM(ROH,Submaps)
  
  if(missing(outfile))
    outfile <- paste("roh_",family.id,"_",individual.id,"_",unit,".png",sep="")
  else{
    outfile <- paste(outfile,".png",sep="") 
  }
  
  plot.segments.id(byROHfile=TRUE, fileOrSubmaps=ROH, individual.id=individual.id, unit=unit, regions=myreg, main=paste("ROHs of ",family.id,"_",individual.id,sep=""), build=build)
}

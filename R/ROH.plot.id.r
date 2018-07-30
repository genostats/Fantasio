#' Creation of plot using ROH file
#' 
#' This function is use to plot HBD segments for a specific infdividual using an ROH file
#' 
#' @param ROHfile an ROH file
#' @param submaps a list.submaps object
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
#' @param regions a matrix containing the value to ve highlighted in the plot
#' @param outfile the name of the plot
#' @param family_id the individual id of the individual wanted
#' @param individual_id the family id of the individual wanted
#' @param save_file whether you want or not to save the plot
#' @param build the value of the build to use to plot chromosome in the plot value accepted are 35, 36, 37, 38 (default is 37)
#' 
#' @details Use this function when you want to use an ROH file to plot your HBD segments for a specific individual.
#' @details Two unit are accepted : "Bases" or "cM".
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end
#'  
#' @return This function returns a plot of the HBD segments for a specific individual
#' 
#' @export
ROH.plot.id <- function(ROHfile, submaps, unit="cM", regions, outfile, family_id, individual_id, save_file=F, build=37)
{
  ROH <- read.table(ROHfile, header=TRUE)
  ROH <- subset(ROH,ROH$FID==family_id & ROH$IID==individual_id)
  
  if(nrow(ROH) == 0)
    return("No information found for this particular individual, please enter a right family_id and individual_id, and make sure individua_id is a character string")
  
  if (missing(regions))
    myreg <- NULL
  else{
    myreg <- regions
  }
  
  if (unit=="cM")
    ROH <- add_cM(ROH,submaps)
  
  if(missing(outfile))
    outfile <- paste("roh_",family_id,"_",individual_id,"_",unit,".png",sep="")
  else{
    outfile <- paste(outfile,".png",sep="") 
  }
  
  if(save_file)
  {
    png(filename = outfile, width = 1000, height = 1000,pointsize=24)
    plot_ROH_IID(ROH,unit, myreg, main=paste("ROHs of ",family_id,"_",individual_id,sep=""), build=build)
    dev.off()
  }
  
  plot_ROH_IID(ROH, unit, myreg, main=paste("ROHs of ",family_id,"_",individual_id,sep=""), build=build)
  
}

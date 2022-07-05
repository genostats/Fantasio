#' Creation of plot using ROH file
#' 
#' This function plots HBDsegments for a specific individual using an ROH file
#' 
#' @param ROH a data frame
#' @param Submaps a atlas object
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
#' @param regions a matrix containing the value to ve highlighted in the plot
#' @param outfile the name of the plot
#' @param id the individual id of the individual wanted
#' @param famid the family id of the individual wanted
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
#' @keywords internal 
plotROHSegmentsId <- function(Submaps, ROH, unit="cM", regions, outfile, famid, id, build)
{
  ROH <- subset(ROH,ROH$FID==famid & ROH$IID==id)
  
  if(nrow(ROH) == 0)
    return("No information found for this individual, please check famid and id values")
  
  if (missing(regions))
    myreg <- NULL
  else{
    #myreg <- regions[regions$chr == chr,]
    myreg <- regions
  }
  
  if (unit=="cM")
    ROH <- addCM(ROH,Submaps)
  
  if(missing(outfile))
    outfile <- paste("roh_",famid,"_",id,"_",unit,".png",sep="")
  else{
    outfile <- paste(outfile,".png",sep="") 
  }
  
  plotSegmentsId(byROHfile=TRUE, fileOrSubmaps=ROH, unit=unit, regions=myreg, main=paste("ROHs of ",famid,":",id,sep=""), build=build)
}

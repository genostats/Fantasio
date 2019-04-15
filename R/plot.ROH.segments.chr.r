#' Creation of plot using ROH file
#' 
#' This function is used to plot HBDsegments for a chromosome using an ROH file
#' 
#' @param ROH a data frame
#' @param submaps a atlas object
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
#' @param chr the chromosome number from which to plot ROH
#' @param outfile the name of the plot
#' @param listid a vector containing the family id follow by the individual id 
#' @param regions a matrix containing the value to ve highlighted in the plot
#' @param build the value of the build to use to plot chromosome in the plot value accepted are 35, 36, 37, 38 (default is 37)
#' 
#' @details Use this function when you want to use an ROH file to plot your HBDsegments for a specific chromosome.
#' @details Two unit are accepted : "Bases" or "cM".
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end 
#' 
#' @return This function returns a plot of the HBDsegments for a specific chromosome
#' @keywords internal 
plot.ROH.segments.chr <- function(ROH, submaps, unit = "cM", chr, outfile, listid, regions, build)
{
  ROH$IID <- as.character(ROH$IID)
  
  if(missing(listid)) 
    listid <- as.character(unique(ROH$IID))[1:10]

  
  if(missing(regions)) 
    myreg <- NULL 
  else {
    myreg <- regions[regions$chr == chr,]
  }
  
  
  
  ROH <- subset(ROH, ROH$CHR==chr)
  
  
  if(unit=="cM") 
    ROH <- add_cM(ROH,submaps)
  
  
  if(missing(outfile))
    outfile <- paste("roh_chr_",chr,".png",sep="")
  else{ 
    outfile <- paste(outfile,".png",sep="") 
  }

  
  plot.segments.chr(byROHfile=TRUE, fileOrSubmaps=ROH,unit= unit,chr= chr,list_id = listid, regions = myreg, build=build)
}

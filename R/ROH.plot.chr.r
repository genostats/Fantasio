#' Creation of plot using ROH file
#' 
#' This function is use to plot HBD segments for a chromosome using an ROH file
#' 
#' @param ROHfile an ROH file
#' @param submaps a list.submaps object
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM")
#' @param chr the chromosome number from which to plot ROH
#' @param outfile the name of the plot
#' @param listid a vector containing the family id follow by the individual id 
#' @param regions a matrix containing the value to ve highlighted in the plot
#' 
#' @details Use this function when you want to use an ROH file to plot your HBD segments for a specific chromosome.
#' @details Two unit are accepted : "Bases" or "cM".
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details -the chromosome number 
#' @details -start 
#' @details -end 
#' 
#' @return This function returns a plot of the HBD segments for a specific chromosome
#' 
#' @export
ROH.plot.chr <- function(ROHfile, submaps, unit = "cM", chr, outfile, listid, regions)
{
  ROH <- read.table(ROHfile, header=TRUE)
  
  ROH$IID <- as.character(ROH$IID)
  if(missing(listid)) 
    list_id <- as.character(unique(ROH$IID))[1:10]
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

  
  plot_ROH_CHR(ROH = ROH, unit = unit, chr = chr, list_id = list_id, regions = myreg)
}

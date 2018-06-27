#' Creation of plot using ROH file
#' 
#' This function is use to plot HBD segments for a chromosome using an ROH file
#' 
#' @param ROHfile an ROH file
#' @param submaps a list.submaps object
#' @param unit the unit used to plot 
#' @param chr the chromosome from which plot will be made
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
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListByHotspots(bedMatrix)
#' submaps <- makeSubmapsByHotspots(bedMatrix, 10, segmentList)  
#' ROH.plot.chr(yourROHfile, submaps)
#' @export
ROH.plot.chr <- function(ROHfile, submaps, unit = "cM", chr, outfile, listid, regions)
{
  ROH <- read.table(ROHfile,h=T)
  
  ROH$IID <- as.character(ROH$IID)
  #for (i in 1:nrow(ROH))
  #  ROH$IID[i] <- paste(ROH$FID[i],ROH$IID[i],sep="_") 
  
  if(missing(listid)) 
    list_id <- as.character(unique(ROH$IID))[1:10]
  #else{ 
  #  table_id <- listid
  #  list_id <- as.character(paste(table_id[,1],table_id[,2],sep="_") )
  #}
  
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

  
  #if(save_file)
  #{
  #  png(filename = outfile, width = 2100, height = 1000,pointsize=24)
  #  #par(mar = c(4.1, 10.1, 4.1, 2.1)); 
  #  plot_ROH_CHR(ROH,unit,chr,list_id,myreg,start=mystart)
  #  dev.off()
  #}
  
  
  plot_ROH_CHR(ROH = ROH, unit = unit, chr = chr, list_id = list_id, regions = myreg)
}

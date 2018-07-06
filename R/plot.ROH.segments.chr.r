plot.ROH.segments.chr <- function(ROHfile, submaps, unit = "cM", chr, outfile, listid, regions)
{
  ROH <- read.table(ROHfile,h=T)
  
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

  
  plot.segments.chr(byROHfile=TRUE, fileOrSubmaps=ROH,unit= unit,chr= chr,list_id = listid, regions = myreg)
}
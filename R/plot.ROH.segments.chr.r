plot.ROH.segments.chr <- function(ROHfile, submaps, unit = "cM", chr, outfile, listid, regions)
{
  ROH <- read.table(ROHfile,h=T)
  
  ROH$IID <- as.character(ROH$IID)
  #for (i in 1:nrow(ROH))
  #  ROH$IID[i] <- paste(ROH$FID[i],ROH$IID[i],sep="_") 
  
  if(missing(listid)) 
    listid <- as.character(unique(ROH$IID))[1:10]
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
  
  
  plot.segments.chr(byROHfile=TRUE, fileOrSubmaps=ROH,unit= unit,chr= chr,list_id = listid, regions = myreg)
}
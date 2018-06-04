ROH.plot.id <- function(ROHfile, submaps, unit="cM", regions, outfile, family_id, individual_id, save_file=F  )
{
  ROH <- read.table(ROHfile,h=T)
  ROH <- subset(ROH,ROH$FID==family_id & ROH$IID==individual_id)
  
  if(nrow(ROH) == 0)
    return("No information found for this particular individual, please enter a right family_id and individual_id, and make sure individua_id is a character string")
  
  if (missing(regions))
    myreg <- NULL
  else{
    myreg <- regions[regions$chr == chr,]
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
    plot_ROH_IID(ROH,unit, myreg, main=paste("ROHs of ",family_id,"_",individual_id,sep=""))
    dev.off()
  }
  
  plot_ROH_IID(ROH, unit, myreg, main=paste("ROHs of ",family_id,"_",individual_id,sep=""))
  
}

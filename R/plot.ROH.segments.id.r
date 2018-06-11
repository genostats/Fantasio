##################################################################################
#This function plots of the HBD probabilities from an ROH file                   #
#                                                                                #
#!!! Submaps : the list of object                                                #                                       
#!!! ROHfile : an ROH file                                                       #
#!!! unit : cM or Bases                                                          #
#!!! regions : a region to be emphasize in the plot                              #
#!!! outfile: (optional) a name for the plot                                     #
#!!! family.id : the family id                                                   #
#!!! individual.id = the individual id                                           #
#                                                                                #
#*** return a plot                                                               #
##################################################################################

plot.ROH.segments.id <- function(Submaps, ROHfile, unit="cM", regions, outfile, family.id, individual.id)
{
  ROH <- read.table(ROHfile,h=T)
  ROH <- subset(ROH,ROH$FID==family.id & ROH$IID==individual.id)
  
  if(nrow(ROH) == 0)
    return("No information found for this particular individual, please enter a right family.id and individual.id, and make sure individua_id is a character string")
  
  if (missing(regions))
    myreg <- NULL
  else{
    myreg <- regions[regions$chr == chr,]
  }
  
  if (unit=="cM")
    ROH <- add_cM(ROH,Submaps)
  
  if(missing(outfile))
    outfile <- paste("roh_",family.id,"_",individual.id,"_",unit,".png",sep="")
  else{
    outfile <- paste(outfile,".png",sep="") 
  }
  
  #if(save_file)
  #{
  #  png(filename = outfile, width = 1000, height = 1000,pointsize=24)
  #  plot_ROH_IID(ROH,unit, myreg, main=paste("ROHs of ",family.id,"_",individual.id,sep=""))
  #  dev.off()
  #}
  
  plot.segments.id(byROHfile=TRUE, fileOrSubmaps=ROH, individual.id=individual.id, unit=unit, regions=myreg, main=paste("ROHs of ",family.id,"_",individual.id,sep=""))
}

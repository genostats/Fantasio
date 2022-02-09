##################################################################################
#This function is uses to create a plot of the HBD for a chromosome and for all  #
#the individual on the submap                                                    #
#                                                                                #
#!!! Submaps : the list of object                                                #                                       
#!!! unit : bases or cM                                                          #
#!!! chr : the number of the chromosome wanted                                   #
#!!! list.ids : (optional) a list of indivial to plot with                       #
#!!! regions  : a regions to be emphasize                                        #
#!!! outfile : the name of the plot                                              #
#                                                                                #
#*** return a new submap object                                                  #
##################################################################################

plotHBDSegmentsChr <- function(Submaps, unit, chr, list.ids, regions, outfile, build)
{
  HBD.recap <- Submaps@HBD_recap
  HBDsegments <- Submaps@HBDsegments
  
  if(missing(list.ids)) 
  {
    # individuals_name <- rownames(HBD.recap)#get the name of the individual
    # individuals_name <- strsplit(individuals_name, "_")
    # individuals_name <- sapply(individuals_name, function(i) match(i[2], Submaps@bedmatrix@ped$id))
    # individuals_name <- paste(Submaps@bedmatrix@ped$famid[individuals_name],Submaps@bedmatrix@ped$id[individuals_name], sep = "_")
    # list.ids <- individuals_name
    list.ids <- rownames(HBD.recap)
  }
    
  
  if(missing(regions)) 
    myreg <- NULL
  else { 
    myreg <- regions[regions$chr == chr,]
  }

  #Get the lines for the wanted chromosome in the HBDsegments dataframe
  HBDsegments_rbind <- do.call(rbind, HBDsegments) #bniding lines
  
  HBD <- subset(HBDsegments_rbind, HBDsegments_rbind$chromosome==chr)#only the wanted lines
  
  HBD$individual <- as.character(HBD$individual) #otherwise factors level in the vector
  HBD$family <- as.character(HBD$family)
  
  #name the file
  if (missing(outfile) )
    outfile <- paste("HBD_chr_",chr,"_",unit,".png",sep="") 
  else 
    outfile <- paste(outfile,".png",sep="") 
  
  plotSegmentsChr(fileOrSubmaps=HBD,unit= unit,chr= chr,list_id = list.ids,regions = myreg, build=build)
}

##################################################################################
#This function plots the HBD probabilities by segments                           #
#                                                                                #
#!!! Submaps : the list of object                                                #                                       
#!!! unit : cM or Bases                                                          #
#!!! regions : a region to be emphasize in the plot                              #
#!!! outfile: (optional) a name for the plot                                     #
#!!! family.id : the family id                                                   #
#!!! individual.id = the individual id                                           #
#                                                                                #
#*** return a plot                                                               #
##################################################################################

plot.HBDsegments.id <- function(Submaps, unit= "cM", individual.id, family.id, regions, outfile, build)
{
  if(!is.character(individual.id))
    return("Need individual id as character")
  if(!is.character(family.id))
    return("Need family id as character")
  
  HBD.recap <- Submaps@HBD_recap
  HBDsegments <- Submaps@HBDsegments

  # individuals_name <- rownames(HBD.recap)#get the name of the individual
  # individuals_name <- strsplit(individuals_name, "_")
  # individuals_name <- sapply(individuals_name, function(i) match(i[2], Submaps@bedmatrix@ped$id))
  # family_id <- Submaps@bedmatrix@ped$famid[individuals_name]
  # individuals_name <- Submaps@bedmatrix@ped$id[individuals_name]
  
  # id   <- which(individuals_name == individual.id)

  # if(length(id) == 0)
  #  stop("No individual found")
  
  HBDsegments_rbind <- do.call(rbind, HBDsegments) #binding lines 
  
  HBD <- subset(HBDsegments_rbind, HBDsegments_rbind$individual==individual.id & HBDsegments_rbind$family==family.id)

  if(nrow(HBD) == 0)
    stop("No individual found")
  
  #regions options
  if (missing(regions)) 
    myreg <- NULL
  else
    myreg <- regions
  
  #name the file
  if (missing(outfile)) 
    outfile <- paste("HBD_", individual.id,"_",unit,".png",sep="")
  else {
    outfile <- paste(outfile,".png",sep="") 
  }

  plot.segments.id(fileOrSubmaps=HBD, unit = unit, regions = myreg, main=paste("HBDsegments of", unique.ids(family.id, individual.id)), build=build)
}


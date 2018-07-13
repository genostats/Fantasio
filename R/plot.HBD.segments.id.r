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

plot.HBD.segments.id <- function(Submaps, unit= "cM", individual.id, family.id, regions, outfile)
{
  if(!is.character(individual.id))
    return("Need individual id as character")
  if(!is.character(family.id))
    return("Need family id as character")
  
  HBD.recap <- Submaps@HBD_recap
  HBD.segments <- Submaps@HBD_segments

  individuals_name <- rownames(HBD.recap)#get the name of the individual
  individuals_name <- strsplit(individuals_name, "_")
  individuals_name <- sapply(individuals_name, function(i) match(i[2], Submaps@bedmatrix@ped$id))
  family_id <- Submaps@bedmatrix@ped$famid[individuals_name]
  individuals_name <- Submaps@bedmatrix@ped$id[individuals_name]
  
  id   <- which(individuals_name == individual.id)
  if(length(id) == 0)
    stop("No individual founded, please make sure that typo is correct.")
  
  HBD_segments_rbind <- do.call(rbind, HBD.segments) #coller toutes les lignes entres elles
  
  HBD <- subset(HBD_segments_rbind, HBD_segments_rbind$individual==individual.id & HBD_segments_rbind$family==family.id)
  if(nrow(HBD) == 0)
    stop("No individual founded, please make sure that typo is correct for individual.id and family.id arguments and in character.")
  
  #traitement de l'options regions
  if (missing(regions)) 
    myreg <- NULL
  else
    myreg <- regions
  
  #donner un nom au fichier creer
  if (missing(outfile)) 
    outfile <- paste("HBD_", individuals_name[id],"_",unit,".png",sep="")
  else {
    outfile <- paste(outfile,".png",sep="") 
  }

  plot.segments.id(fileOrSubmaps=HBD, individual.id=id, unit = unit, regions = myreg, main=paste("HBD segments of ",family.id, "_", individual.id, sep = ""))
}


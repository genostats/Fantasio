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

plot.HBD.segments.chr <- function(Submaps, unit, chr, list.ids, regions, outfile)
{
  HBD.recap <- Submaps@HBD_recap
  HBD.segments <- Submaps@HBD_segments
  
  if(missing(list.ids)) 
  {
    individuals_name <- rownames(HBD.recap)#get the name of the individual
    individuals_name <- strsplit(individuals_name, "_")
    individuals_name <- sapply(individuals_name, function(i) match(i[2], Submaps@bedmatrix@ped$id))
    individuals_name <- paste(Submaps@bedmatrix@ped$famid[individuals_name],"_",Submaps@bedmatrix@ped$id[individuals_name])
    list.ids <- individuals_name
  }
    
  
  if(missing(regions)) 
    myreg <- NULL
  else { 
    myreg <- regions[regions$chr == chr,]
  }
  #Partie pour recuperer les lignes pour le chromosome voulue dans le dataframe HBD.segments
  HBD_segments_rbind <- do.call(rbind, HBD.segments) #coller toutes les lignes entres elles
  
  HBD <- subset(HBD_segments_rbind, HBD_segments_rbind$chromosome==chr)#recuperer seulement les lignes pour le chromosome voulut
  
  HBD$individual <- as.character(HBD$individual) #si cette etape n'est pas faite les nom des individus sont des factors
  
  #donner un nom au fichier creer
  if (missing(outfile) )
    outfile <- paste("HBD_chr_",chr,"_",unit,".png",sep="") 
  else 
    outfile <- paste(outfile,".png",sep="") 
  
  plot.segments.chr(fileOrSubmaps=HBD,unit= unit,chr= chr,list_id = list.ids,regions = myreg)
}

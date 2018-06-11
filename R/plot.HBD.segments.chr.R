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
  
  #nom <- rownames(HBD.recap)	
  
  #si aucun vecteur list.ids n'est fournit
  if(missing(list.ids)) 
  {
    individuals_name <- rownames(HBD.recap)#get the name of the individual
    individuals_name <- strsplit(individuals_name, "_")
    individuals_name <- sapply(individuals_name, function(i) match(i, Submaps@bedmatrix@ped$id))
    individuals_name <- individuals_name[!is.na(individuals_name)]
    individuals_name <- Submaps@bedmatrix@ped$id[individuals_name]
    
    list.ids <- individuals_name
  }
    
  
  #Traitement de l'options regions
  if(missing(regions)) 
    myreg <- NULL
  else { # data frame avec une colonne chr et 
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
  
  
  #creation du fichier png et appel de la fonction de plot
  #if(save_img)
  #{
  #  png(filename =  outfile, width = 2100, height = 1000,pointsize=24)
  #  #par(mar = c(4.1, 10.1, 4.1, 2.1)); 
  #  plot_HBD_CHR(HBD_segments = HBD,unit= unit,chr= chr,list_id = list.ids,regions = myreg)
  #  dev.off()
  #}
  #par(mar = c(4.1, 10.1, 4.1, 2.1)); 
  plot.segments.chr(fileOrSubmaps=HBD,unit= unit,chr= chr,list_id = list.ids,regions = myreg)
}

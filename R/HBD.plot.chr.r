##################################################################################
#This function allows to plot graphs about the HBD probabilities for a chromosome#
#and every individual                                                            #
#!!! Submaps : the list of objects                                               #                                       
#!!! ROHfile : (optionnal) if you want to create the plot using an ROH file      #
#!!! unit    : bases or cM  3                                                    #
#!!! chr     : the number of the chromosome you want                             #
#!!! list.ids: (optional) a list of individual you want                          #
#!!! regions : a region you want to be emphasize                                 #
#!!! outfile : the name of the plot                                              #
#                                                                                #
#*** return a plot                                                               #
##################################################################################

#HBD.recap : dataframe avec les moyennes des probabilites HBD pour tous les individus et tout les marqueurs apparu au moins une fois
#HBD.segments : list de dataframe pour chaque individus avec les segments contenant au moins 5 marqueurs consecutifs au dessus ou egale a un threshold
#chr : le numero du chromosome
#listid : la liste des individus a prendre en compte 

HBD.plot.chr <- function(Submaps, ROHfile, unit="cM", chr, list.ids, regions, outfile)
{
  if(!missing(Submaps) & !missing(ROHfile))
  {
    plot.ROH.segments.chr(ROHfile = ROHfile, submaps = Submaps, unit = unit, chr = chr, outfile=outfile, listid=list.ids, regions=regions)
  }else{
    if(!missing(Submaps))
      plot.HBD.segments.chr(Submaps=Submaps, unit=unit, chr=chr, list.ids=list.ids, regions=regions, outfile=outfile)
  }
}
  
  
  
  #HBD.recap <- Submaps@HBD_recap
  #HBD.segments <- Submaps@HBD_segments
  #
  #nom <- rownames(HBD.recap)	
  ##si aucun vecteur list.ids n'est fournit
  #if(missing(list.ids)) 
  #  list.ids <- as.character(unique(nom)) 
 #
  #  #Traitement de l'options regions
  #if(missing(regions)) 
  #  myreg <- NULL
  #else { # data frame avec une colonne chr et 
  #  myreg <- regions[regions$chr == chr,]
  #}
  #
  ##Partie pour recuperer les lignes pour le chromosome voulue dans le dataframe HBD.segments
  #HBD_segments_rbind <- do.call(rbind, HBD.segments) #coller toutes les lignes entres elles
  #
  #HBD <- subset(HBD_segments_rbind, HBD_segments_rbind$chromosome==chr)#recuperer seulement les lignes pour le chromosome voulut
  #
  #HBD$individual <- as.character(HBD$individual) #si cette etape n'est pas faite les nom des individus sont des factors
  #
  ##donner un nom au fichier creer
  #if (missing(outfile) )
  #  outfile <- paste("HBD_chr_",chr,"_",unit,".png",sep="") 
  #else 
  #  outfile <- paste(outfile,".png",sep="") 
  #
  #
  ##creation du fichier png et appel de la fonction de plot
  ##if(save_img)
  ##{
  ##  png(filename =  outfile, width = 2100, height = 1000,pointsize=24)
  ##  #par(mar = c(4.1, 10.1, 4.1, 2.1)); 
  ##  plot_HBD_CHR(HBD_segments = HBD,unit= unit,chr= chr,list_id = list.ids,regions = myreg)
  ##  dev.off()
  ##}
  ##par(mar = c(4.1, 10.1, 4.1, 2.1)); 
  #plot_HBD_CHR(HBD_segments = HBD,unit= unit,chr= chr,list_id = list.ids,regions = myreg)
  
  



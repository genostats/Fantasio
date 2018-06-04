##################################################################################
#This function create a plot for an individuals                                  #
#                                                                                #
#!!! Submaps : a list object                                                     #                                       
#!!! ROHfile : (optional) a ROHfile                                              #
#!!! unit : cM ou Bases                                                          #
#!!! individual.id : the individual id                                           #
#!!! individual.id : the individual family id                                    #
#!!! regions : (optional) a regions you want to be emphasize                     #
#!!! outfille : (optional) the name of the plot                                  #
#                                                                                #
#*** return a plot                                                               #
##################################################################################

HBD.plot.id <- function(Submaps, ROHfile, unit= "cM", individual.id, family.id, regions, outfile)
{
  if(!missing(Submaps) & !missing(ROHfile))
  {
    plot.ROH.segments.id(Submaps=Submaps, ROHfile, unit, individual.id=individual.id, family.id=family.id, regions, outfile=outfile)
  }else{
    if(!missing(Submaps))
      plot.HBD.segments.id(Submaps = Submaps, individual.id=individual.id, family.id=family.id, unit=unit, regions, outfile=outfile)
  }
}

  
  #if(!is.character(individual.id))
  #  return("Need individual.id as character")
  #
  #HBD.recap <- Submaps@HBD_recap
  #HBD.segments <- Submaps@HBD_segments
  #
  ##recuperer les id de l'individus necessaire pour la fonction plot
  #name <- rownames(HBD.recap)
  #id   <- which(name == individual.id)
  #
  ##traitement de l'options regions
  #if (missing(regions)) 
  #  myreg <- NULL
  #else
  #  myreg <- regions
  #
  ##donner un nom au fichier creer
  #if (missing(outfile)) 
  #  outfile <- paste("HBD_", name[id],"_",unit,".png",sep="")
  #else {
  #  outfile <- paste(outfile,".png",sep="") 
  #}
#
  ##creation d'un fichier png et plot
  ##if(save_img)
  ##{
  ##  png(filename = outfile, width = 1000, height = 1000,pointsize=24)
  ##  plot_HBD_IID(HBD.segments, id, unit = unit, regions = myreg, main=paste("HBD segments of ",name[id]))
  ##  dev.off()
  ##}
  #
  #plot_HBD_IID(HBD.segments, id, unit = unit, regions = myreg, main=paste("HBD segments of ",name[id]))

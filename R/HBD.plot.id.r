#HBD.recap : dataframe avec les moyennes des probabilites HBD pour tous les individus et tout les marqueurs apparu au moins une fois
#HBD.segments : list de dataframe pour chaque individus avec les segments contenant au moins 5 marqueurs consecutifs au dessus ou egale a un threshold
#individual.number : le numero de l'individus sortis dans HBD.recap ou HBD.segments (meme indidividus) dont on veut le plot

HBD.plot.id <- function(HBD.recap, HBD.segments, individual.number, distance= "cM", regions="empty", outfile="empty")
{
  #recuperer les id de l'individus necessaire pour la fonction plot
  name <- rownames(HBD.recap)
  
  #traitement de l'options regions
  if (regions=="empty") {myreg=NULL;} else {myreg=read.table(regions,h=F);}
  
  #donner un nom au fichier creer
  if (outfile=="empty") {outfile=paste("HBD_", name[individual.number],"_",distance,".png",sep="")} else { outfile=paste(outfile,".png",sep="") }

  #creation d'un fichier png et plot
  png(filename = outfile, width = 1000, height = 1000,pointsize=24)
  plot_HBD_IID(HBD.segments, individual.number, distance = distance, regions = myreg, main=paste("HBD segments of ",name[individual.number]))
  dev.off()
}



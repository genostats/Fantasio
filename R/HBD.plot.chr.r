#HBD.recap : dataframe avec les moyennes des probabilites HBD pour tous les individus et tout les marqueurs apparu au moins une fois
#HBD.segments : list de dataframe pour chaque individus avec les segments contenant au moins 5 marqueurs consecutifs au dessus ou egale a un threshold
#chr : le numero du chromosome
#listid : la liste des individus a prendre en compte 

HBD.plot.chr <- function(HBD.recap, HBD.segments, distance="cM", chr, listid="empty", regions="empty", outfile="empty")
{
  
  nom <- rownames(HBD.recap)	
  #si aucun fichier lis_id n'est fournit
  if (listid=="empty") 
  {
    list_id=as.character(unique(nom)) 
  }
  
  #si un fichier list_id est fournit
  else
  { 
    table_id=read.table(listid,h=F)
    list_id=as.character(paste(table_id[,1],table_id[,2],sep="_") )
  }
  
  #Traitement de l'options regions
  if (regions=="empty") {myreg=NULL;}
  else {myreg=read.table(regions,h=F); myreg=myreg[myreg[,1]==chr,]}
  
  #Partie pour recuperer les lignes pour le chromosome voulue dans le dataframe HBD.segments
  HBD_segments_rbind <- do.call(rbind, HBD.segments) #coller toutes les lignes entres elles
  
  HBD=subset(HBD_segments_rbind, HBD_segments_rbind$chromosome==chr)#recuperer seulement les lignes pour le chromosome voulut
  
  HBD$individual <- as.character(HBD$individual) #si cette etape n'est pas faite les nom des individus sont des factors
  
  #donner un nom au fichier creer
  if (outfile=="empty") {outfile=paste("HBD_chr_",chr,"_",distance,".png",sep="")} else { outfile=paste(outfile,".png",sep="") }
  
  #special cases for chromosomes when genetic maps starts after 10 Mb
  mystart=""
  if (distance=="cM"){
    if( (chr==13) || (chr==14) || (chr==15)){
      mystart=-20
    } else if ((chr==21) || (chr==22)) {
      mystart=-15
    }
  } 
  
  #creation du fichier png et appel de la fonction de plot
  png(filename = outfile, width = 2100, height = 1000,pointsize=24)
  par(mar = c(4.1, 10.1, 4.1, 2.1)); 
  plot_HBD_CHR(HBD_segments = HBD,distance= distance,chr= chr,list_id = list_id,regions = myreg,start=mystart)
  dev.off()

}


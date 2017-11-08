plot_HBD_CHR <- function (HBD_segments, distance = "cM", chr, list_id, regions, start="", end="", color2="green4")
{
  
  options(echo=FALSE);
  suppressMessages(library(quantsmooth))
  
  color  <- function (val) {
    if (val==1) {col="skyblue3"}
    else if (val==2) {col="tomato"}
    else {col="grey"}
    col
  }
  
  #choisir dans quelle unite le plot sera fait
  
                               #numero de colonne debut segment                  numero de colonne fin du segment dans le dataframe HBD_segments
  if (distance=="cM")         {pos1=which(colnames(HBD_segments)=="start_dist"); pos2=which(colnames(HBD_segments)=="end_dist"); myxlab="Position (cM)"; coeff=1; if (length(pos1)==0) {stop("no genetic distances in HBD file")}}
  else if (distance=="bases") {pos1=which(colnames(HBD_segments)=="start_pos"); pos2=which(colnames(HBD_segments)=="end_pos"); myxlab="Position (Mb)"; coeff=1000000}
  else {stop("length option accepts cM and bases only")}
  
  
  if (start=="") {start=0}
  if (end=="")   {end=lengthChromosome(chr,distance)}
  
  div=length(list_id)/10
  
  #creer un plot vide 
  plot(c(start,end),c(-0.5,length(list_id))/div+0.5,type="n",xaxt="n",yaxt="n",ylab="",xlab=myxlab,main=paste("HBD segments on chromosome ",chr,sep=""), cex.main=1.5, cex.lab=1.5,font.lab=2)
  #	axis(1, at = (end-abs(start))/2 ,labels = myxlab,col.ticks="white",font.axis=2, cex.axis=1.5 )
  
  
  #traitement de l'option regions
  if (is.null(regions)==F){
    if (nrow(regions)>0) {         #dessiner les regions
      for (i in 1:nrow(regions)) {polygon(regions[i,c(2,3,3,2)],c(rep(0.5,2),rep(length(list_id)+0.5,2))/div,col=color2,border=color2,lwd=2)}
    }
  }
  
  #dessiner le chromosome
  paintCytobands(chr,units=distance,pos=c(0,0),orientation="h",legend = FALSE)
  for (j in 1:length(list_id)){
    
    axis(2,at=(length(list_id)-j+1)/div,list_id[j],col.ticks=0,las=2,font=2,cex.axis=1.25)
    
    toplot=HBD_segments[HBD_segments$individual==list_id[j],]
    
    if (nrow(toplot) >0 ) {
      for (k in 1:nrow(toplot)) {
        polygon(c(toplot[k,pos1],toplot[k,pos2],toplot[k,pos2],toplot[k,pos1]),c(rep(length(list_id)-j+1.25,2),rep(length(list_id)-j+0.75,2))/div,col=color(toplot$status[k]),lwd=2)
      }
    }
  }
  
  cpt=0
  
  while (cpt*20*coeff < end) {
    axis(1,at=cpt*20*coeff,cpt*20,cex.axis=1.25)
    cpt=cpt+1
  }
  
}

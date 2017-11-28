plot_HBD_IID <- function (HBD_segments, individual_number, distance = "cM", regions="empty",color2="green4", main="") 
{
  
  options(echo=FALSE);
  suppressMessages(library(quantsmooth))
  
  color  <- function (val) {
    if(is.na(val)){col="black"}
    else if (val==1) {col="skyblue3"}
    else if (val==2) {col="tomato"}
    else {col="grey"}
    col
  }
  
  #choisir dans quelle unite le plot sera fait
  
  if (distance=="cM")         {ecart=25; larg=2; pos1=which(colnames(HBD_segments[[individual_number]])=="start_dist"); pos2=which(colnames(HBD_segments[[individual_number]])=="end_dist"); if (length(pos1)==0) {stop("no genetic distances for this individual")}}
  else if (distance=="bases") {ecart=25000000; larg=2; pos1=which(colnames(HBD_segments[[individual_number]])=="start_pos"); pos2=which(colnames(HBD_segments[[individual_number]])=="end_pos");}
  else {stop("length option accepts cM and bases only")}
  
  #creer un plot vide
  plot(c(larg*0.5,larg*11.5),c(0,lengthChromosome(1,distance)+lengthChromosome(12,distance)+1.25*ecart),type="n",xaxt="n",yaxt="n",xlab="",ylab="",main=main,cex.main=1.5)
  
  #on fait deux boucle pour pouvoir avoir deux rangees de chromosome, premiere de 1 a 11 et deuxieme de 12 a 22
  
  #premiere boucle pour creer les 11 premiers chromosomes
  
  for (i in 1:11) {
    #dessiner le chromosome
    paintCytobands(i,units=distance,pos=c(larg*i,lengthChromosome(12,distance)+lengthChromosome(i,distance)+1.25*ecart),orientation="v",legend = FALSE)
    
    text(larg*i,lengthChromosome(12,distance)+0.75*ecart,i)
    
    #recuperer les lignes pour un individus specifique et un chromosome specifique
    seg_chr=HBD_segments[[individual_number]][HBD_segments[[individual_number]]$chromosome==i,] 
    
    if (dim(seg_chr)[1]>0)
    {
        for (k in 1:dim(seg_chr)[1]) 
        {
          xx=c(rep(larg*(i+0.00),2),rep(larg*(i+0.5),2))
          
          yy=c(lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr[k,pos1],seg_chr[k,pos2])),lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr[k,pos2],seg_chr[k,pos1])))
          
          polygon(xx,yy,col=color(seg_chr$status[k]),border=color(seg_chr$status[k]))
          #lines(rep(larg*(i+0.5),2),lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr$deb[k],seg_chr$fin[k])),lwd=5, col=col_freq(seg_chr$freq[k]))
        }
    }
    
    #traitement de l'option regions
    reg_chr=regions[regions[,1]==i,]
    
    if (is.null(reg_chr)==F)
    {
      if(dim(reg_chr)[1]>0) 
      {
        for (k in 1:dim(reg_chr)[1]) 
        {
          xx=c(rep(larg*(i+0.00),2),rep(larg*(i+0.5),2))
            yy=c(lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr[k,2],reg_chr[k,3])),lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr[k,3],reg_chr[k,2])))
            polygon(xx,yy,border=color2,lwd=2)
            #lines(rep(larg*(i+0.5),2),lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr$deb[k],reg_chr$fin[k])),lwd=5, col=col_freq(reg_chr$freq[k]))
        }
      }
    }
  }
  
  
  #deuxieme boucle pour creer les 11 derniers chromosomes
  
  for (i in 12:22) 
  {
    #dessiner le chromosome
    paintCytobands(i,units=distance,pos=c(larg*(i-11),lengthChromosome(i,distance)+ecart/4),orientation="v",legend = FALSE)
    
    text(larg*(i-11),-0.25*ecart,i)
    
    #recuperer les lignes pour un individus et un chromosome specifique
    seg_chr=HBD_segments[[individual_number]][HBD_segments[[individual_number]]$chromosome==i,] 
    
    if (dim(seg_chr)[1]>0) 
    {
      for (k in 1:dim(seg_chr)[1]) 
      {
        xx=c(rep(larg*(i-11+0.00),2),rep(larg*(i-11+0.5),2))
        
        yy=c(0.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr[k,pos1],seg_chr[k,pos2])),0.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr[k,pos2],seg_chr[k,pos1])))
        
        polygon(xx,yy,col=color(seg_chr$status[k]),border=color(seg_chr$status[k]))
        #lines(rep(larg*(i-11+0.5),2),ecart/4+abs(lengthChromosome(i,distance)-c(seg_chr$deb[k],seg_chr$fin[k])),lwd=5, col=col_freq(seg_chr$freq[k]))
      }
      
    }
    
    #traitement de l'options regions
    reg_chr=regions[regions[,1]==i,]
    
    if (is.null(reg_chr)==F)
    {
      if(dim(reg_chr)[1]>0) 
      {
        for (k in 1:dim(reg_chr)[1]) 
        {
          xx=c(rep(larg*(i-11+0.00),2),rep(larg*(i-11+0.5),2))
          
          yy=c(0.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr[k,2],reg_chr[k,3])),0.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr[k,3],reg_chr[k,2])))
          
          polygon(xx,yy,border=color2,lwd=2)
          #lines(rep(larg*(i-11+0.5),2),ecart/4+abs(lengthChromosome(i,distance)-c(reg_chr$deb[k],reg_chr$fin[k])),lwd=5, col=col_freq(reg_chr$freq[k]))
        }
      }
    }
    
  }
  
}

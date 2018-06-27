#' Plot of the HFLOD 
#' 
#' This fonction plot the HFLOD score for a chromosome
#' 
#' @param submaps a list.submaps object
#' @param unit the unit used to plot 
#' @param chr the chromosome number from which to plot HFLOD score
#' @param regions a matrix containing the value to ve highlighted in the plot
#' @param color2 the color of the regions highlighted
#' @param nbSNP_MA
#' 
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details - the chromosome number 
#' @details - start 
#' @details - end
#' 
#' 
#' @return This function returns a manhattan plot of all the HFLOD score over all the chromosome
#' 
#' 
#' @examples  
#' bedMatrix <- read.bed.matrix("yourFile")
#' segmentList <- createSegmentsListByHotspots(bedMatrix)
#' submaps <- makeSubmapsByHotspots(bedMatrix, 10, segmentList)  
#' HFLOD.manhattan.plot(h)
#' @export
HFLOD.plot.chr <- function(submaps, unit = "cM", chr, regions, color2="green4", nbSNP_MA = 50) 
{
  if(submaps@bySegments)
  {
    HFLOD <- submaps@HFLOD
    pos <- submaps@atlas[[1]]@map$dist
    chromosome <- submaps@bedmatrix@snps$chr[submaps@atlas[[1]]@submap]
  }
  #HFLOD=read.table(paste(folder,"/HFLOD.temp.txt",sep=""),h=T)
  #pos <- sapply(submaps@atlas, function(hh) hh@submap)
  #pos <- unique(unlist(pos))
  else{
    HFLOD <- submaps@HFLOD
    pos <- as.data.frame(table(unlist(sapply(submaps@atlas, function(hh) hh@submap))), stringsAsFactors=FALSE)
    pos <- as.numeric(pos$Var1)
    pos <- sort(pos)
    chromosome <- submaps@bedmatrix@snps$chr[pos]
    pos <- submaps@bedmatrix@snps$dist[pos]
  }
  
  

  
  #if (regions!="empty") {myreg=read.table(regions,h=F); color2="green4"}
  if(missing(regions)) 
    myreg <- NULL
  else { 
    myreg <- regions
    color2="green4"
    myreg$start = regions$start/1e6
    myref$end   = regions$end/1e6
  }
  
  
  #if (distance=="cM"){
  #  pos    <- 3
  #  myxlab <- "Position (cM)"
  #} else {
  #  pos <- 4
  #  HFLOD[,pos]=HFLOD[,pos]/1000000;                               
  #  myxlab="Position (Mb)"                               
  #  if (regions!="empty") {myreg[,2:3]=myreg[,2:3]/1000000}
  #}
  
  if(unit == "cM"){
    myxlab <- "Position (cM)"
    coeff  <- 1
  }else{
    myxlab <- "Position (Mb)"
    coeff  <- 1e6
    pos    <- pos/1e6
  }
  
  #1)graphs per chromosome
  newout   <- NULL
  axis_mp  <- NULL
  chr_pos  <- 5 
  myreg_mp <- NULL
  
  
  #for (c in 1:unique(submaps@atlas[[1]]@map$chr)){
    #h@atlas[[1]]@HFLOD[h@atlas[[1]]@map$chr == 1, 1]
    
    
    #toplot <- HFLOD[HFLOD$CHR==c,]
    toplot_HFLOD <- HFLOD[chromosome == chr, 1]
    toplot_MA    <- HFLOD[chromosome == chr, 2]
    toplot_pos   <- pos[chromosome == chr]
    
    ymax <- max(3.3,max(toplot_HFLOD))
    
    #png(file=paste(folder,"/HFLOD.",distance,".",c,".png",sep=""), width = 2400, height = 1000,pointsize=24)
    
    #plot(toplot[,pos],toplot$HFLOD,pch=16,ylim=c(0,ymax),xlab=myxlab,ylab="HFLOD",cex.lab=1.4,cex.axis=1.5,main=paste("HFLOD (chromosome ",c,")",sep=""),cex.main=1.5)
    
    plot(toplot_pos,
         toplot_HFLOD,
         pch  = 16,
         ylim = c(0,ymax),
         xlab = myxlab,
         ylab = "HFLOD",
         cex.lab  = 1.4,
         cex.axis = 1.5,
         main     = paste("HFLOD (chromosome ",chr,")",sep=""),
         cex.main = 1.5)


    #if (regions!="empty") {
    #  myreg_chr=subset(myreg,myreg[,1]==c)
    #  if (nrow(myreg_chr)>0) {
    #    for (i in 1:nrow(myreg_chr)) {
    #      polygon(myreg_chr[i,c(2,3,3,2)],c(rep(-1,2),rep(ymax+1,2)),col=color2,border=color2,lwd=2)
    #      myreg_mp=rbind(myreg_mp,max(c(0,axis_mp))+10++myreg_chr[i,2:3])
    #    }
    #    points(toplot[,pos],toplot$HFLOD,pch=16)
    #  }
    #}

    if(!(missing(regions))){
      myreg_chr <-  myreg[which(myreg$CHR == chr),]
      if(nrow(myreg_chr) > 0){
        for(i in 1:nrow(myreg_chr)){
          #polygon(myreg_chr[i,c(2,3,3,2)],c(rep(-1,2),rep(ymax+1,2)),col=color2,border=color2,lwd=2)
          polygon(x = myreg_chr[i,c("start","end","end","start")]/coeff, 
                  y = c(rep(-1,2),rep(ymax+1,2)), 
                  col    = color2,
                  border = color2,
                  lwd    = 2)
          
          myreg_mp = rbind(myreg_mp,max(c(0,axis_mp))+10+myreg_chr$start[i]+myreg_chr$end[i])
        }
        points(toplot_pos, toplot_HFLOD, pch=16)
      }
    }
    for (i in 1:3) 
      abline(h=i,col="grey",lwd=1,lty=2)

    abline(h=3.3,col="grey",lwd=2)
    
    #lines(toplot[,pos],rollmean(toplot$HFLOD,as.numeric(nbSNP_MA),fill="extend"), col="red",lwd=5)
    lines(toplot_pos, zoo::rollmean(toplot_HFLOD, as.numeric(nbSNP_MA), fill = "extend"), col="red", lwd=2)
    
    axis_mp <- c(axis_mp, max(c(0,axis_mp))+10+toplot_pos)
    chr_pos <- c(chr_pos, max(c(0,axis_mp))+5)
    
    #newout=rbind(newout,cbind(toplot,round(rollmean(toplot$HFLOD,as.numeric(nbSNP_MA),fill="extend"),3)))
    #newout <- rbind(newout,cbind(toplot,round(rollmean(toplot$HFLOD,as.numeric(nbSNP_MA),fill="extend"),3)))
    
  #}
  #colnames(newout)=c("CHR","RS","POS_cM","POS_bp","HFLOD","ALPHA","MA_HFLOD")
  #if (distance=="bases") {newout[,pos]=newout[,pos]*1000000;}
  #HFLOD=newout
  #write.table(HFLOD,file=paste(folder,"/HFLOD.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
  
  
  
  
  
  # #2)Manhattan plot
  # png(file=paste(folder,"/HFLOD.",distance,".png",sep=""), width = 2400, height = 1000,pointsize=24)
  # plot(axis_mp,HFLOD$HFLOD,pch=16,ylim=c(0,ymax),xlab="",ylab="HFLOD",cex.lab=1.4,cex.axis=1.5,col=HFLOD$CHR,xaxt="n",cex=0.75)
  # for(i in 1:22) {abline(v=chr_pos[i],col="grey"); axis(1,at=mean(chr_pos[i:(i+1)]),i,col.ticks=0,cex.axis=1.5)}
  # abline(v=chr_pos[23],col="grey");
  # abline(h=3,col="grey")
  # dev.off()
  
  # #2)Manhattan plot smooth
  # png(file=paste(folder,"/HFLOD.MA.",distance,".png",sep=""), width = 2400, height = 1000,pointsize=24)
  # plot(axis_mp,HFLOD$MA_HFLOD,pch=16,ylim=c(0,ymax),xlab="",ylab="HFLOD with moving average",cex.lab=1.4,cex.axis=1.5,col=HFLOD$CHR,xaxt="n",cex=0.75)
  # for(i in 1:22) {abline(v=chr_pos[i],col="grey"); axis(1,at=mean(chr_pos[i:(i+1)]),i,col.ticks=0,cex.axis=1.5)}
  # abline(v=chr_pos[23],col="grey");
  # abline(h=3,col="grey")
  # dev.off()
  
  
}

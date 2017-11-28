add_cM <- function (ROH, map) {
	
	#add 3 colones "POS1_cM", "POS2_cM" and "cM"
	ROH=cbind(ROH,rep(0,nrow(ROH)),rep(0,nrow(ROH),nrow(ROH)),rep(0,nrow(ROH),nrow(ROH)))
	colnames(ROH)[ncol(ROH)-2]="POS1_cM"; colnames(ROH)[ncol(ROH)-1]="POS2_cM"; colnames(ROH)[ncol(ROH)]="cM"
	
	#read the map with cM positions
	if (map=="") {stop("no map file specified\n")}
	else         {map=read.table(map)}	
	
	#add cM information
	for (i in 1:nrow(ROH)){
		ROH$POS1_cM[i]=subset(map$V3,as.character(map$V2)==as.character(ROH$SNP1[i]))
		ROH$POS2_cM[i]=subset(map$V3,as.character(map$V2)==as.character(ROH$SNP2[i]))
		ROH$cM[i]     =ROH$POS2_cM[i]-ROH$POS1_cM[i]
	}
		
	ROH

}

color  <- function (val) {
  if(is.na(val)){col="black"}
	else if (val==1) {col="skyblue3"}
	else if (val==2) {col="tomato"}
	else {col="grey"}
	col
}

plot_ROH_IID <- function (ROH, distance, regions,color2="green4", main="") {

	if (distance=="cM")         {ecart=25; larg=2; pos1=which(colnames(ROH)=="POS1_cM"); pos2=which(colnames(ROH)=="POS2_cM"); if (length(pos1)==0) {stop("no genetic distances in ROH file")}}
	else if (distance=="bases") {ecart=25000000; larg=2; pos1=which(colnames(ROH)=="POS1")   ; pos2=which(colnames(ROH)=="POS2");}
	else {stop("length option accepts cM and bases only")}
	
  #plot vide
	plot(c(larg*0.5,larg*11.5),c(0,lengthChromosome(1,distance)+lengthChromosome(12,distance)+1.25*ecart),type="n",xaxt="n",yaxt="n",xlab="",ylab="",main=main,cex.main=1.5)
	
	#11 premiers chromosomes
	for (i in 1:11) {
		paintCytobands(i,units=distance,pos=c(larg*i,lengthChromosome(12,distance)+lengthChromosome(i,distance)+1.25*ecart),orientation="v",legend = FALSE)
		text(larg*i,lengthChromosome(12,distance)+0.75*ecart,i)
		seg_chr=ROH[ROH$CHR==i,]
		if (dim(seg_chr)[1]>0) {
		for (k in 1:dim(seg_chr)[1]) {
			xx=c(rep(larg*(i+0.00),2),rep(larg*(i+0.5),2))
			yy=c(lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr[k,pos1],seg_chr[k,pos2])),lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr[k,pos2],seg_chr[k,pos1])))
			polygon(xx,yy,col=color(seg_chr$PHE[k]),border=color(seg_chr$PHE[k]))
			#lines(rep(larg*(i+0.5),2),lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr$deb[k],seg_chr$fin[k])),lwd=5, col=col_freq(seg_chr$freq[k]))
		}
		}
		reg_chr=regions[regions[,1]==i,]
		if (is.null(reg_chr)==F){
		if (dim(reg_chr)[1]>0) {
			for (k in 1:dim(reg_chr)[1]) {
				xx=c(rep(larg*(i+0.00),2),rep(larg*(i+0.5),2))
				yy=c(lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr[k,2],reg_chr[k,3])),lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr[k,3],reg_chr[k,2])))
				polygon(xx,yy,border=color2,lwd=2)
				#lines(rep(larg*(i+0.5),2),lengthChromosome(12,distance)+1.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr$deb[k],reg_chr$fin[k])),lwd=5, col=col_freq(reg_chr$freq[k]))
			}
		}
		}
	}
	
	#11 derniers chromosomes
	for (i in 12:22) {
		paintCytobands(i,units=distance,pos=c(larg*(i-11),lengthChromosome(i,distance)+ecart/4),orientation="v",legend = FALSE)
		text(larg*(i-11),-0.25*ecart,i)
		seg_chr=ROH[ROH$CHR==i,]
		if (dim(seg_chr)[1]>0) {
		for (k in 1:dim(seg_chr)[1]) {
			xx=c(rep(larg*(i-11+0.00),2),rep(larg*(i-11+0.5),2))
			yy=c(0.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr[k,pos1],seg_chr[k,pos2])),0.25*ecart+abs(lengthChromosome(i,distance)-c(seg_chr[k,pos2],seg_chr[k,pos1])))
			polygon(xx,yy,col=color(seg_chr$PHE[k]),border=color(seg_chr$PHE[k]))
			#lines(rep(larg*(i-11+0.5),2),ecart/4+abs(lengthChromosome(i,distance)-c(seg_chr$deb[k],seg_chr$fin[k])),lwd=5, col=col_freq(seg_chr$freq[k]))
		}
		}
		reg_chr=regions[regions[,1]==i,]
		if (is.null(reg_chr)==F){
		if (dim(reg_chr)[1]>0) {
			for (k in 1:dim(reg_chr)[1]) {
				xx=c(rep(larg*(i-11+0.00),2),rep(larg*(i-11+0.5),2))
				yy=c(0.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr[k,2],reg_chr[k,3])),0.25*ecart+abs(lengthChromosome(i,distance)-c(reg_chr[k,3],reg_chr[k,2])))
				polygon(xx,yy,border=color2,lwd=2)
				#lines(rep(larg*(i-11+0.5),2),ecart/4+abs(lengthChromosome(i,distance)-c(reg_chr$deb[k],reg_chr$fin[k])),lwd=5, col=col_freq(reg_chr$freq[k]))
			}
		}
		}
		
	}

}

ROH.plot.id <- function(ROHfile, map, distance="cM", regions="empty", outfile="empty", fid, iid, save_file=F  )
{
  ROH=read.table(ROHfile,h=T)
  ROH=subset(ROH,ROH$FID==fid & ROH$IID==iid)
  if (regions=="empty") { myreg=NULL } else {myreg=read.table(regions,h=F)}
  if (distance=="cM") {ROH=add_cM(ROH,map)}
  if (outfile=="empty") {outfile=paste("roh_",fid,"_",iid,"_",distance,".png",sep="")} else { outfile=paste(outfile,".png",sep="") }
  if(save_file)
  {
    png(filename = outfile, width = 1000, height = 1000,pointsize=24)
    plot_ROH_IID(ROH,distance, myreg, main=paste("ROHs of ",fid,"_",iid,sep=""))
    dev.off()
  }
  plot_ROH_IID(ROH, distance, myreg, main=paste("ROHs of ",fid,"_",iid,sep=""))
  
}


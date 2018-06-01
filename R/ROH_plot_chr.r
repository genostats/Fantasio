add_cM <- function (ROH, map) {
	
	#add 3 colones "POS1_cM", "POS2_cM" and "cM"
	ROH=cbind(ROH,rep(0,nrow(ROH)),rep(0,nrow(ROH),nrow(ROH)),rep(0,nrow(ROH),nrow(ROH)))
	colnames(ROH)[ncol(ROH)-2]="POS1_cM"; colnames(ROH)[ncol(ROH)-1]="POS2_cM"; colnames(ROH)[ncol(ROH)]="cM"
	
	#read the map with cM positions
	if (map=="") {stop("no map file specified\n")}
	else         {map=read.table(map)}	
	
	#add cM information
	if(nrow(ROH)>0){
		for (i in 1:nrow(ROH)){
			ROH$POS1_cM[i]=subset(map$V3,as.character(map$V2)==as.character(ROH$SNP1[i]))
			ROH$POS2_cM[i]=subset(map$V3,as.character(map$V2)==as.character(ROH$SNP2[i]))
			ROH$cM[i]     =ROH$POS2_cM[i]-ROH$POS1_cM[i]
		}
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
####################################################################################################################
plot_ROH_CHR <- function (ROH, distance, chr, list_id, regions, start="", end="", color2="green4") {

	if (distance=="cM")         {pos1=which(colnames(ROH)=="POS1_cM"); pos2=which(colnames(ROH)=="POS2_cM"); myxlab="Position (cM)"; coeff=1; if (length(pos1)==0) {stop("no genetic distances in ROH file")}}
	else if (distance=="bases") {pos1=which(colnames(ROH)=="POS1")   ; pos2=which(colnames(ROH)=="POS2"   ); myxlab="Position (Mb)"; coeff=1000000;}
	else {stop("length option accepts cM and bases only")}

	if (start=="") {start=0}
	if (end=="")   {end=lengthChromosome(chr,distance)}
	
	div=length(list_id)/10;
	#plot vide
	plot(c(start,end),c(-0.5,length(list_id))/div+0.5,type="n",xaxt="n",yaxt="n",ylab="",xlab=myxlab,main=paste("ROHs on chromosome ",chr,sep=""), cex.main=1.5, cex.lab=1.5,font.lab=2)
#	axis(1, at = (end-start)/2 ,labels = paste("Position (",distance,")",sep=""),col.ticks="white",font.axis=2, cex.axis=1.5 )
	
	#option regions
	if (is.null(regions)==F)
	{ 
	  if (nrow(regions)>0) 
	    {
		    for (i in 1:nrow(regions)) {polygon(regions[i,c(2,3,3,2)],c(rep(0.5,2),rep(length(list_id)+0.5,2))/div,col=color2,border=color2,lwd=2)}
	    }
	}
	#dessiner le chromosome
	paintCytobands(chr,units=distance,pos=c(0,0),orientation="h",legend = FALSE)
	for (j in 1:length(list_id)){
		axis(2,at=(length(list_id)-j+1)/div,list_id[j],col.ticks=0,las=2,font=2,cex.axis=1.25)
		toplot=ROH[ROH$IID==list_id[j],]
		if (nrow(toplot) >0 ) {
			for (k in 1:nrow(toplot)) {
				polygon(c(toplot[k,pos1],toplot[k,pos2],toplot[k,pos2],toplot[k,pos1]),c(rep(length(list_id)-j+1.25,2),rep(length(list_id)-j+0.75,2))/div,col=color(toplot$PHE[k]),lwd=2)
			}
		}
	}
	cpt=0
	while (cpt*20*coeff < end) {
		axis(1,at=cpt*20*coeff,cpt*20,cex.axis=1.25)
		cpt=cpt+1
	}
}

#####################################################################################################################
ROH.plot.chr <- function(ROHfile, map, chr, outfile, listid="empty", regions="empty", distance="cM", save_file=F)
{
  ROH=read.table(ROHfile,h=T)
  
  ROH$IID=as.character(ROH$IID); for (i in 1:nrow(ROH)){ ROH$IID[i]=paste(ROH$FID[i],ROH$IID[i],sep="_") }	
  if (listid=="empty") { list_id=as.character(unique(ROH$IID)) 
  } else { 
  	table_id=read.table(listid,h=F)
  	list_id=as.character(paste(table_id[,1],table_id[,2],sep="_") )
  }
  if (regions=="empty") { myreg=NULL } else {myreg=read.table(regions,h=F); myreg=myreg[myreg[,1]==chr,]}
  ROH=subset(ROH,ROH$CHR==chr)
  if (distance=="cM") {ROH=add_cM(ROH,map)}
  if (outfile=="empty") {outfile=paste("roh_chr_",chr,".png",sep="")} else { outfile=paste(outfile,".png",sep="") }
  
  #special cases for chromosomes when genetic maps starts after 10 Mb
  mystart=""
  if (distance=="cM"){
  	if( (chr==13) || (chr==14) || (chr==15)){
   		mystart=-20
   	} else if ((chr==21) || (chr==22)) {
   		mystart=-15
   	}
  }
  
  if(save_file)
  {
    png(filename = outfile, width = 2100, height = 1000,pointsize=24)
    par(mar = c(4.1, 10.1, 4.1, 2.1)); 
    plot_ROH_CHR(ROH,distance,chr,list_id,myreg,start=mystart)
    dev.off()
  }
  plot_ROH_CHR(ROH, distance,chr,list_id,myreg,start=mystart)
}


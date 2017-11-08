#ROH=cbind(ROH,rep(0,nrow(ROH)),rep(0,nrow(ROH),nrow(ROH)),rep(0,nrow(ROH),nrow(ROH)))
#colnames(ROH)[ncol(ROH)-2]="POS1_cM"; colnames(ROH)[ncol(ROH)-1]="POS2_cM"; colnames(ROH)[ncol(ROH)]="cM"

#paintCytobands(chr,units=distance,pos=c(0,0),orientation="h",legend = FALSE)
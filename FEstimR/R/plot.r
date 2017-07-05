
n[which(n$INBRED_Fsuite == 1 & n$INBRED_FEstimR == 1),]  # meme inbred
n[which(n$INBRED_Fsuite == 1 & n$INBRED_FEstimR == 0),]  #inbred pour Fsuite mais pas pour FestimR
n[which(n$INBRED_Fsuite == 0 & n$INBRED_FEstimR == 1),]  #inbre pour FestimR mais pas pour Fsuite
n[which(n$INBRED_Fsuite == 0 & n$INBRED_FEstimR == 0),]  #aucun inbred



#par(xpd=TRUE, mar=c(5, 4, 4, 12.5))
#plot(n$Fsuite_Fmean, n$FestimR_Fmean, xlab="Fsuite", ylab="Festim", main="Comparaison F_Mean", col=c("black", "red", "green", "cyan")[n$coltest])
#namevector <- c("Both INBRED", "FSuite INBRED / FEstimR NOT INBRED", "FEstimR INBRED / FSuite NOT INBRED", "NO INBRED")
#legend("topright",inset=c(-0.65,0), legend= namevector, col=c("black", "red", "green", "cyan"), pch=1, cex=0.7)
#par(xpd=FALSE)
#abline(0,1,col="red")


par(xpd=TRUE, mar=c(5, 4, 4, 12.5))
plot(FSuite$F_MEDIAN, FEstimR$F_MEDIAN, xlab="Fsuite", ylab="FestimR", main="Comparaison F_Mean", col=c("black", "red", "green", "cyan")[n$coltest])
namevector <- c("Both INBRED", "FSuite INBRED / FEstimR NOT INBRED", "FEstimR INBRED / FSuite NOT INBRED", "NO INBRED")
legend("topright",inset=c(-0.65,0), legend= namevector, col=c("black", "red", "green", "cyan"), pch=1, cex=0.7)
par(xpd=FALSE)
abline(0,1,col="red")

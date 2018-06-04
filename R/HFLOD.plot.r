##plot for a known chromosom
##
#HFLOD.plot.chr.msat.matrix <- function(msatmatrix, chr) 
#{
#  #find the distance
#  w <- which(msatmatrix@map$chr == chr)
#  d <- msatmatrix@map$distance[w]
#  n <- length(d)
#  r <- msatmatrix@HFLOD[w, 1]
#  
#  #plotting
#  ymax <- max(3.3, max(msatmatrix@HFLOD[w,], na.rm = TRUE)) 
#  xmax <- max(d, na.rm = TRUE)
#  ymin <- 0
#  #d=distance for each chr i ; r=HFLOD for each individual j & chr i
#  plot(d, r, type="l", pch=16, xlim = c(0,xmax), ylim = c(ymin, ymax), 
#      xlab = "Position (cM)", ylab="HFLOD", main= paste("HFLOD (chromosome ",chr,")", sep=""))    
#  
#  for ( m in 1:ymax)
#    abline(h=m, col="grey", lwd=1.5, lty=2)
#  
#  abline(h=3.3, col="grey", lwd=2)
#}

###############################################
##saving plot for all the chromosom on png format
#
#HFLOD.plot.chrs.png <- function(x)
#{ 
#
#  dir.create("Plots")#create de Plots folder
#  for ( i in unique(x@map$chr)) # for each chr we make a plot 
#  {
#    png(file = paste("Plots/HFLOD.","cM.",i,".png",sep = ""), width = 2400, height = 1000, pointsize=24)
#    par(cex = 1.45)
#    HFLOD.plot.chr(x,i)
#    dev.off()
#  }  
#}
#
###############################################
#manhattan plot 
#HFLOD.manhattan.plot.msat.matrix <- function(x)
#{
#  #take the last position of the chr i
#  res <- NULL
#  res_pos <- 5
#  
#  for ( i in unique(x@map$chr))
#  {
#    d <- x@map$distance[x@map$chr == i]
#    res <- c(res, max(c(0,res, na.rm=TRUE))+10+d)
#    res_pos <- c(res_pos, max(c(0,res))+5)
#  }
#  
#  # other way to do 
#  
#  #tapply( x@map$distance, x@map$chr, function(d) tail(d,1) ) -> chr_le
#  #tapply( x@map$distance, x@map$chr, length ) -> chr_n
#  #nb_chr <- length(chr_le)
#  #offset <- c(10, cumsum(chr_le+10)[ -nb_chr ] )
#  #pos <- x@map$distance + unlist(mapply(rep, offset, chr_n))
#  #res_pos <- c(6,cumsum(chr_le +10)+6)
#  
#  ymax <- max(3.3,max(x@HFLOD[,1], na.rm = TRUE))
#  ymin <- 0
#  xmax <- max(res, na.rm=TRUE)
#  mycol <- rep(c("cadetblue2",8),11) #couleur des points
#  plot(res, x@HFLOD[,1],pch=16,xlim=c(0,xmax), ylim = c(ymin,ymax), xlab = "", ylab = "HFLOD", cex.lab = 1.4, cex.axis=1.5, main= paste("HFLOD : Manhattan Plot ", sep=""), col= mycol[x@map$chr], xaxt="n", cex=0.75 )  
#  #pour chaque chr les separer par des lignes
#  for ( j in unique(x@map$chr))
#  {
#    abline(v=res_pos[j],col="grey",lwd=2)
#    axis(1,at=mean(res_pos[j:(j+1)]),j,col.ticks=0,cex.axis=1.5)
#  }
#  for (j in 1:3) 
#  {
#    abline(h=j,col="grey",lwd=1,lty=2)
#  }
#  abline(h=3.3,col="grey",lwd=2)
#}
# 
#  
#  
#  
#
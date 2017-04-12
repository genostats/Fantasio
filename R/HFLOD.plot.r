HFLOD.plot <- function(x)
{ 
  
  ###############################
  #       graph per chr         #
  ###############################
  
  
  dir.create("Plots")#creation du dossier Plots contenant les sorties
  last <- function(x) { return (length(x))}
  
  
  for ( i in unique(x@map$chr)) # pour chaque chr on fait un plot 
  {
    r <- NULL
    #Obtenir les distances
    d <- x@map$distance[x@map$chr == i]
    n <- last(d)
    r <- x@HFLOD[ x@map$chr == i,1]
    
    
    #plotting
    ymax <- max(3, max(x@HFLOD[x@map$chr == i, ], na.rm = TRUE)) 
    xmax <- max(d, na.rm = TRUE)
    ymin <- min(x@HFLOD[x@map$chr == i,], na.rm = TRUE)
    png(file = paste("Plots/HFLOD.","cM.",i,".png",sep = ""), width = 4000, height = 2500, pointsize=24)
    
    #d=distance pour chaque chr i ; r=HFLOD pour chaque Indiv j & chr i
    plot(d, r, type="l", pch=16, xlim = c(0,xmax), ylim = c(ymin, ymax), xlab = "Position (cM)", ylab="HFLOD", cex.lab=1.4, cex.axis=1.5, main= paste("HFLOD (chromosome ",i,")", sep=""), cex.main=1.5)    
    
    for ( m in 1:ymax)
    {
      abline(h=m, col="grey", lwd=1.5, lty=2)
    }
    
    abline(h=3, col="grey", lwd=2)
    
    
    #Obtenir les HFLOD de chaque individu
    #for ( j in 1:x@nrow )
    #{
      #if(x@f[j] == 0) next
      
      #lines(d, r,col= j)
      
    #} 
    legend(x = "topright", legend =  c("Individus :", x@ped$famid, "Famid : ",  x@ped$id), lty =0,cex = 1.5) 
    
    
    
    dev.off()
  }  
  
  ############################
  #      Manhattan plot      #
  ############################
  
  #prendre la derniere position du chr i
  
  res <- NULL
  res_pos <- 5
  
  for ( i in unique(x@map$chr))
  {
    d <- x@map$distance[x@map$chr == i]
    res <- c(res, max(c(0,res, na.rm=TRUE))+10+d)
    res_pos <- c(res_pos, max(c(0,res))+5)
  }
  

  
  ymax <- max(3,max(x@HFLOD[,1], na.rm = TRUE))
  ymin <- min(x@HFLOD, na.rm=TRUE)
  xmax <- max(res, na.rm=TRUE)
  mycol <- rep(c("cadetblue2",8),11) #couleur des points
  
  #white mont plot 
  
  #if(x@f[i] == 0) next
  png(file = paste("Plots/HFLOD.cM.png",sep=""), width = 2400, height = 800, pointsize = 24 )
  plot(res, x@HFLOD[,1],pch=16,xlim=c(0,xmax), ylim = c(ymin,ymax), xlab = "", ylab = "HFLOD", cex.lab = 1.4, cex.axis=1.5, main= paste("HFLOD : White Mont Plot ", sep=""), col= mycol[x@map$chr], xaxt="n", cex=0.75 )  
  
  #pour chaque chr les separer par des lignes
  for ( j in unique(x@map$chr))
  {
    abline(v=res_pos[j],col="grey",lwd=2)
    axis(1,at=mean(res_pos[j:(j+1)]),j,col.ticks=0,cex.axis=1.5)
  }
  for (j in 1:3) 
  {
    abline(h=j,col="grey",lwd=1,lty=2)
  }
  abline(h=3,col="grey",lwd=2)
  dev.off()
  
  
  
  
}
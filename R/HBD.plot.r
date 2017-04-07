HBD.plot <- function(x)
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
      
      
      #plotting
      ymax <- max(3, max(x@FLOD[, x@map$chr == i], na.rm = TRUE)) 
      xmax <- max(d, na.rm = TRUE)
      ymin <- min(x@FLOD[, x@map$chr == i], na.rm = TRUE)
      png(file = paste("Plots/FLOD.","cM.",i,".png",sep = ""), width = 2400, height = 1000, pointsize=24)
          
      #d=distance pour chaque chr i ; r=FLOD pour chaque Indiv j & chr i
      plot(d, r, type="n", pch=16, xlim = c(0,xmax), ylim = c(ymin, ymax), xlab = "Position (cM)", ylab="FLOD", cex.lab=1.4, cex.axis=1.5, main= paste("FLOD (chromosome ",i,")", sep=""), cex.main=1.5)    
      
      for ( m in 1:ymax)
      {
        abline(h=m, col="grey", lwd=1.5, lty=2)
      }
      
      abline(h=3, col="grey", lwd=2)
      
      
      #Obtenir les FLOD de chaque individu
      for ( j in 1:x@nrow )
        
      {
        if(x@f[j] == 0) next
        r <- x@FLOD[j, x@map$chr == i]
        lines(d, r,col= j)
      } 
      
      
                  
      
      
      dev.off()
  }  
    
  ############################
  #      Manhattan plot      #
  ############################
  
  #prendre la derniere position du chr i
  n <- 0
  res <- c(x@map$distance[x@map$chr == 1])
  res_pos <- c()
  
  for ( i in unique(x@map$chr))
  {
    d <- x@map$distance[x@map$chr == i]
    n2 <- last(d)
    
   
    # ajouter a toutes les positions du chr i + 1 la derniere valeurs 
    #de la pos du chr i
   
    d <- x@map$distance[x@map$chr == i+1] + n
    v <- d[last(n2)]
    n <- n + x@map$distance[n2]
    res <- c(res,d)
    res_pos <- c(res_pos, v)
  }
  res_pos <- res_pos[!is.na(res_pos)]
  last_res <- last(res)
  h <- last(res_pos)
  res_pos[h +1] <- res[last_res] # utile pour tracer les lignes
  
  
  ymax <- max(x@FLOD, na.rm = TRUE)
  ymin <- min(x@FLOD, na.rm=TRUE)
  xmax <- max(res, na.rm=TRUE)
  mycol <- rep(c("cadetblue2",8),11) #couleur des points
  
  #Plotting : tous les chr pour un individus 
  png(file = paste("Plots/FLOD.", "cM.","png",sep=""), width = 2400, height = 800, pointsize = 24 )
  for ( i in x@nrow)
  {
    if(x@f[j] == 0) next
    plot(res, x@FLOD[i,], type="p", pch=16, xlim=c(0,xmax), ylim = c(ymin,ymax), xlab = "", ylab = "FLOD", cex.lab = 1.4, cex.axis=1.5, col= mycol[unique(x@map$chr)], xaxt="n", cex=0.75 )  #
    #lines(res,x@FLOD[i,],col=2,lwd=3) # demande des valeurs finie  
    #for(i in unique(x@map$chr)) 
    #{
    #  abline(v=x@map$distance[],col="grey",lwd=2);
    #  axis(1,at=mean(chr_pos[i:(i+1)]),i,col.ticks=0,cex.axis=1.5)
    #}
    #abline(v=chr_pos[23],col="grey",lwd=2);  #
    #for ( i in 1:ymax)
    #{
    #  abline(h=i, col="grey", lwd=1.5, lty=2)
    #}  
    #abline(h=3,col="grey",lwd=2)
    
    
    #pour chaque chr les separer par des lignes
    n <- NULL
    for ( j in unique(x@map$chr))
    {
      d <- x@map$distance[x@map$chr == j]
      n <- last(d)
      d <- x@map$distance[x@map$chr == j+1] + n
      
      axis(1,at=res_pos,1:22,col.ticks=0,cex.axis=1.5)
      abline(v=res_pos,col="grey",lwd=2)
      n <- n + x@map$distance[n2]
    }
     
  }
  
  dev.off()
    
}
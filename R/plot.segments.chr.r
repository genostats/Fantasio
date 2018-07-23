plot.segments.chr <- function(byROHfile=FALSE, fileOrSubmaps, unit = "cM", chr, list_id, regions, color2="green4", build=37)
{
  if(length(list_id) > 20 )
    list_id <- list_id[1:20]
  
  l <- unit.plot.chr(file = fileOrSubmaps, unit, byROHfile = byROHfile) 
  pos1   <- l$pos1
  pos2   <- l$pos2
  myxlab <- l$myxlab
  coeff  <- l$coeff
  
  #special cases for chromosomes when genetic maps starts after 10 Mb
  start <- 0
  if (unit=="cM"){
    if( (chr==13) || (chr==14) || (chr==15)){
      start <- -20
    } else if ((chr==21) || (chr==22)) {
      start <- -15
    }
  } 
  
  end <- lengthChromosome(chr,unit, build)/coeff
  
  
  #creer un plot vide 
  y_max <- length(list_id)+1
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, y_max/2.5, 0, 0)) 
  
  plot(x <- c(start,end), y <- c(0,y_max), 
       type="n", yaxt="n", ylab="", xlab=myxlab, 
       main=paste("HBD segments on chromosome ",chr,sep=""))

  axis(2, at = (1:length(list_id))+0.25, list_id, col.ticks=0, las=2) 
  
  #traitement de l'option regions
  if(!is.null(regions)){
    if (nrow(regions)>0) {         #dessiner les regions
      for (i in 1:nrow(regions)) {
        polygon(x = regions[i,c("start","end","end","start")]/coeff,
                y = c(0.75,0.75,y_max,y_max), col=color2, border=color2, lwd=2)
      }
    }
  }
  
  #dessiner le chromosome
  paintCytobands(chr,units=unit,pos=c(0,0.5), build=build, orientation="h",legend = FALSE, length.out = end)
  
  #dessiner les segments HBD
  tmp <- strsplit(list_id, "_")
  list_id <- sapply(tmp, function(i) i[2])
  
  for (j in 1:length(list_id)){
    if(byROHfile)
    {
      toplot <- fileOrSubmaps[fileOrSubmaps$IID==list_id[j],]
    }else{
      toplot <- fileOrSubmaps[fileOrSubmaps$individual==list_id[j],]
    }
    
    
    if (nrow(toplot) > 0) {
      for (k in 1:nrow(toplot)) {
        polygon( x  = c(toplot[k,pos1],toplot[k,pos2],toplot[k,pos2],toplot[k,pos1])/coeff,
                 y  = c(j,j,j+0.5,j+0.5),
                 col=ifelse(byROHfile, color(toplot$PHE[k]), color(toplot$status[k])),
                 lwd=1)
      }
    }
  }
  
  
}

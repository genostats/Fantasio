plot.segments.id <- function (byROHfile=FALSE, fileOrSubmaps, individual.id, unit = "cM", regions, color2="green4", main) 
{
  #choisir dans quelle unite le plot sera fait
  
  #if (unit=="cM"){
  #  ecart <- 25
  #  larg  <- 2
  #  pos1  <- which(colnames(HBD_segments[[individual_number]])=="start_dist")
  #  pos2  <- which(colnames(HBD_segments[[individual_number]])=="end_dist")
  #  if (length(pos1)==0)
  #    stop("no genetic distances for this individual")
  #}
  #else if (unit=="bases"){
  #  ecart <- 2.5e7
  #  larg  <- 2
  #  pos1  <- which(colnames(HBD_segments[[individual_number]])=="start_pos")
  #  pos2  <- which(colnames(HBD_segments[[individual_number]])=="end_pos")
  #}
  #else{
  #  stop("units parameter accepts cM and bases only")
  #}
  l <- unit.plot.id(file=fileOrSubmaps, unit=unit, byROHfile, individual.id)
  ecart <- l$ecart
  larg  <- l$larg
  pos1  <- l$pos1
  pos2  <- l$pos2
  
  #creer un plot vide
  if(missing(main)) 
    main = NULL
  plot(x = c(larg*0.5,larg*11.5),
       y = c(0,quantsmooth::lengthChromosome(1,unit)+quantsmooth::lengthChromosome(12,unit)+1.25*ecart),
       type="n",
       xaxt="n",
       yaxt="n",
       xlab="",
       ylab="",
       main=main,
       cex.main=1.5)
  
  
  #boucle sur les chr
  for (i in 1:22) {
    if(i < 12) {
      i2 <- i
      offset_y <- quantsmooth::lengthChromosome(12, unit) + ecart
    } else {
      i2 <- i - 11
      offset_y <- 0
    }

    #dessiner le chromosome
    quantsmooth::paintCytobands(i,units=unit,pos=c(larg*i2, offset_y + quantsmooth::lengthChromosome(i,unit) + ecart/4), orientation="v", legend = FALSE)
    
    text(larg*i2, offset_y - 0.25*ecart,i)
    
    #traitement de l'options regions
    if (!is.null(regions)) {
      regions_chr <- regions[ regions$chr == i, ]
      if(nrow(regions_chr)>0) 
      {
        for (k in 1:nrow(regions_chr)) 
        {
          xx <- c(rep(larg*i2,2),rep(larg*(i2+0.5),2))
                  
                  yy <- c(0.25*ecart+abs(quantsmooth::lengthChromosome(i,unit)-c(regions_chr$start[k],regions_chr$end[k])),
                          0.25*ecart+abs(quantsmooth::lengthChromosome(i,unit)-c(regions_chr$end[k],regions_chr$start[k]))) 
                  
                  polygon(x = xx,
                          y = yy + offset_y,
                          col = color2,
                          border = color2,
                          lwd = 1.5)
                  #lines(rep(larg*(i2+0.5),2),ecart/4+abs(lengthChromosome(i,distance)-c(reg_chr$deb[k],reg_chr$fin[k])),lwd=5, col=col_freq(reg_chr$freq[k]))
        }
      }
    }
    
    #recuperer les lignes pour un individus et un chromosome specifique
    if(byROHfile)
    {
      seg_chr <- fileOrSubmaps[fileOrSubmaps$CHR==i,]
    }else{
      seg_chr <- fileOrSubmaps[fileOrSubmaps$chromosome==i,]
    }
    
    if (nrow(seg_chr) > 0) 
    {
      for (k in 1:nrow(seg_chr)) 
      {
        xx <- c(rep(larg*i2,2),rep(larg*(i2+0.5),2))
        
        yy <- c(0.25*ecart+abs(quantsmooth::lengthChromosome(i,unit)-c(seg_chr[k,pos1],seg_chr[k,pos2])),
                0.25*ecart+abs(quantsmooth::lengthChromosome(i,unit)-c(seg_chr[k,pos2],seg_chr[k,pos1])))
        
        polygon(x      = xx,
                y      = yy + offset_y,
                col    = ifelse(byROHfile,color(seg_chr$PHE[k]),color(seg_chr$status[k])),
                border = ifelse(byROHfile,color(seg_chr$PHE[k]),color(seg_chr$status[k])), 
                lwd    = 1)
        #lines(rep(larg*(i2+0.5),2),ecart/4+abs(lengthChromosome(i,distance)-c(seg_chr$deb[k],seg_chr$fin[k])),lwd=5, col=col_freq(seg_chr$freq[k]))
      }
    }
  }
}
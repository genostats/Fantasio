#' Plot of the HFLOD 
#' 
#' This fonction plot the HFLOD score for a chromosome
#' 
#' @param submaps a list.submaps object
#' @param unit the unit used to plot, two options are allowed "Bases", "cM" (default is "CM") 
#' @param chr the chromosome number from which to plot HFLOD score
#' @param regions a matrix containing the value to ve highlighted in the plot
#' @param color2 the color of the regions highlighted (default is "green4")
#' @param nbSNP_MA number of SNP for the moving average (default is 50)
#' 
#' @details If you use the regions options make sure to pass a matrix containing one line per region to be highlighted with in each line : 
#' @details - the chromosome number 
#' @details - start 
#' @details - end
#' 
#' @seealso set.HFLOD
#' 
#' @return This function returns a manhattan plot of all the HFLOD score over all the chromosome
#' 
#' 
#' @examples  
#' ##install.packages("HGDP.CEPH", repos="https://genostats.github.io/R/") ## make this only one time
#' require(Fantasio)
#' require(HGDP.CEPH)
#' filepath <-system.file("extdata", "hgdp_ceph.bed", package="HGDP.CEPH")
#' x <- read.bed.matrix(filepath)
#' x <- set.stats(x)
#' x.me <- select.inds(x, population == "Bedouin")
#' x.me@ped$pheno <- rep(2,48) #The package analyzes only individualw with a status of 2
#' submaps <- Fantasio(x.me, "Hotspots", n=5)
#' HFLOD.manhattan.plot(submaps)
#' @export
HFLOD.plot.chr <- function(submaps, unit = "cM", chr, regions, color2="green4", nbSNP_MA = 50) 
{
  if(class(submaps@bedmatrix)[1] != "bed.matrix")
    stop("Need a bed.matrix to eat")
  
  if(class(submaps@atlas[[1]])[1] != "snps.matrix" & class(submaps@atlas[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.")
  
  if(is.null(submaps@HFLOD))
    stop("HFLOD slots in the object is empty, cannot plot")
  
  if(submaps@bySegments && class(submaps@atlas[[1]])[1] == "snps.matrix")
    stop("Cannot plot by segments for snps.matrix object")
  
  HFLOD <- submaps@HFLOD
  #to get mean position when working by segments
  if(unit == "cM")
    pos <- HFLOD$pos_cM
  else
    pos <- HFLOD$pos_Bp
  
  chromosome <- HFLOD$CHR
  
  
  
  if(missing(regions)) 
    myreg <- NULL
  else { 
    myreg <- regions
    color2="green4"
    myreg$start = regions$start/1e6
    myref$end   = regions$end/1e6
  }
  
  
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
  
  toplot_HFLOD <- HFLOD$HFLOD[HFLOD$CHR == chr]
  toplot_MA    <- HFLOD$HFLOD[HFLOD$CHR == chr]
  if(unit == "cM")
    toplot_pos   <- HFLOD$pos_cM[HFLOD$CHR == chr]
  else
    toplot_pos   <- HFLOD$pos_Bp[HFLOD$CHR == chr]
  
  ymax <- max(3.3,max(toplot_HFLOD))
  
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
  
  
  if(!(missing(regions))){
    myreg_chr <-  myreg[which(myreg$CHR == chr),]
    if(nrow(myreg_chr) > 0){
      for(i in 1:nrow(myreg_chr)){
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
  
  lines(toplot_pos, zoo::rollmean(toplot_HFLOD, as.numeric(nbSNP_MA), fill = "extend"), col="red", lwd=2)
  
  axis_mp <- c(axis_mp, max(c(0,axis_mp))+10+toplot_pos)
  chr_pos <- c(chr_pos, max(c(0,axis_mp))+5)
}

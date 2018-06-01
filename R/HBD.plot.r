#x = bedmatrix
#HBD.plot.chr <- function(HBD_recap, bedmatrix,n.consecutive.marker = 5, threshold = 0.5, chr) 
#{
#  HBD_segments <- HBD.segments(HBD_recap, n.consecutive.marker = 5, threshold = 0.5)#list
#  
#  color  <- function (val) {
#    if (val==1) {col="skyblue3"}
#    else if (val==2) {col="tomato"}
#    else {col="grey"}
#    col
#  }
#  
#  #trouver les marqueurs du chromosome
#  
#  lll <- match(colnames(HBD_recap), bedmatrix@snps$id) #trouver les indices des marqueurs de HBD_prob dans la bedmatrix
#  chromosome <- bedmatrix@snps$chr[lll] #trouver a quelle chromosome appartient le marqueur
#  marqueur_chr <- which(chromosome == chr)#sort les indices du vecteur chromosome
#  val <- HBD_recap[marqueur_chr]
#  dist <- bedmatrix@snps$dist[marqueur_chr]
#  
#  #plotting
#  ymax <- max(1, max(val, na.rm = TRUE)) 
#  xmax <- max(dist, na.rm = TRUE)
#  ymin <- 0
#  #d=distance for each chr i ; r=HFLOD for each individual HBD_recap & chr i
#  plot(dist, val, pch=16, xlim = c(0,xmax), ylim = c(ymin, ymax), 
#       xlab = "Position (cM)", ylab="HBD", main= paste("HBD (chromosome ",chr,")", sep=""))    
#  
#  for ( m in seq(ymin, ymax, by=0.1))
#    abline(h=m, col="grey", lwd=1.5, lty=2)
#  
#}
#
#####################################################
#HBD.plot.ind <- function(HBD_recap, bedmatrix, ind) 
#{
#  #trouver les marqueurs du chromosome
#  
#  lll <- match(colnames(HBD_recap), bedmatrix@snps$id) #trouver les indices des marqueurs de HBD_recap dans la bedmatrix
#  dist <- bedmatrix@snps$dist[lll]
#  
#  #plotting
#  ymax <- max(1, max(HBD_recap, na.rm = TRUE)) 
#  xmax <- max(dist, na.rm = TRUE)
#  ymin <- 0
#  #d=distance for each chr i ; r=HFLOD for each individual HBD_recap & chr i
#  plot(dist, HBD_recap[ind,], pch=16, xlim = c(0,xmax), ylim = c(ymin, ymax), 
#       xlab = "Position (cM)", ylab="HBD", main= paste("HBD (individus ",ind,")", sep=""))    
#  
#  for ( m in seq(ymin, ymax, by=0.1))
#    abline(h=m, col="grey", lwd=1.5, lty=2)
#  
#}#
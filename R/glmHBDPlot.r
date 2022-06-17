#' QQ-Plot & Manhattan Plots for glm on HBD prob or FLOD
#' 
#' @param x an atlas object
#' @param expl_var the explanatory variable 'FLOD' or 'HBD_prob'
#' @param n.cores number of cores for parallelization calculation (default = 1)
#' @param save choose if plot are saved or not (.png format) (default = FALSE)
#' 
#' @seealso \code{\link{glmHBD}}
#' 
#' @export

glmHBDPlot = function ( x, expl_var, save = FALSE) {
	
	unadj 	<- x@logisticRegression$unadj

	adj 		<- x@logisticRegression$adj
	
	adj <- adj [which(adj$p_value != 0), ]
	unadj <- unadj [which(unadj$p_value != 0), ]
	
	
	# QQ Plot
	# unadjusted
	#p1 <- ggQQPlot(unadj$p_value) + ggtitle(paste("QQ-plot GLM with", expl_var, "- unadjusted data"))
	
	# adjusted
	#p2 <- ggQQPlot(adj$p_value) + ggtitle(paste("QQ-plot GLM with", expl_var, "- adjusted data"))
	
	# Manhattan Plot
	# Change colnames to fit with gaston
	colnames(unadj)[colnames(unadj) == 'pos_Bp'] <- 'pos'
	colnames(unadj)[colnames(unadj) == 'p_value'] <- 'p'
	
	colnames(adj)[colnames(adj) == 'pos_Bp'] <- 'pos'
	colnames(adj)[colnames(adj) == 'p_value'] <- 'p'
	
	treshold = -log10(0.05/dim(x@submaps_list[[1]])[2])  #treshold based on the number of markers
	lim <- round(max(-log10(adj$p), -log10(unadj$p), treshold))+1
	
	if (save == FALSE) { # Default, just print plots on different windows

	#opar <- par(mfrow = c(2,1))  
	
	# QQ Plot
	#dev.new(height = 5.1, width = 5)
	gaston::qqplot.pvalues(unadj$p, main = paste("QQ-plot GLM with", expl_var, "- unadjusted data"), ylim = c(0,lim) )
	
	#dev.new(height = 5.1, width = 5, ypos = 650)
	gaston::qqplot.pvalues(adj$p, main = paste("QQ-plot GLM with", expl_var, "- adjusted data"), ylim = c(0,lim))
	
	#dev.new(width = 14.25, height= 5.1, ypos = 0, xpos = 640)
	gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
	abline(h=treshold, col = 'red')
	
	# adjusted
	#dev.new(width = 14.25, height= 5.1, ypos = 650, xpos = 640)
	gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
	abline(h=treshold, col = 'red')
	
	#par(opar)
	
	} else { # save plots in png files 
	
	png(paste('qqplot.GLM.', expl_var,'unadj.png')) 		# U = unadjusted
	gaston::qqplot.pvalues(unadj$p, main = paste("QQ-plot GLM with", expl_var, "- unadjusted data"), ylim = c(0,lim) )
	dev.off()
	
	png(paste('qqplot.GLM.', expl_var, 'adj.png'))		# A = adjusted
	gaston::qqplot.pvalues(adj$p, main = paste("QQ-plot GLM with", expl_var, "- adjusted data"), ylim = c(0,lim))
	dev.off()
	
	png(paste( 'manhattanplot.GLM.', expl_var,'unadj.png'))
	gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
	abline(h=treshold, col = 'red')
	dev.off()	
	
	png(paste( 'manhattanplot.GLM.', expl_var,'adj.png'))
	gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"), ylim = c(0,lim))
	abline(h=treshold, col = 'red')
	dev.off()	
	
	}
			  
		
}



#' QQ-Plot & Manhattan Plots for glm on HBD prob or FLOD
#' 
#' @param x the bedmatrix
#' @param expl_var the explanatory variable 'FLOD' or 'HBD_prob'
#' @param covar_df a dataframe containing covariates
#' @param covar covariates of interest 
#' if missing, all covariates of the dataframe are considered
#' @param n.cores number of cores for parallelization calculation (default = 1)
#' @param save choose if plot are saved or not (default = FALSE)
#' 
#' @seealso glmHBD()
#' 
#' @export

glmHBDPlot = function ( x, expl_var, save = FALSE) {
	
	unadj 	<- x@logisticRegression$unadj

	adj 		<- x@logisticRegression$adj
	
	adj <- adj [which(adj$p_value != 0), ]
	unadj <- unadj [which(unadj$p_value != 0), ]
	
	# QQ Plot
	# unadjusted
	p1 <- ggQQPlot(unadj$p_value) + ggtitle(paste("QQ-plot GLM with", expl_var, "- unadjusted data"))
	
	# adjusted
	p2 <- ggQQPlot(adj$p_value) + ggtitle(paste("QQ-plot GLM with", expl_var, "- adjusted data"))
	
	# Manhattan Plot
	# Change colnames to fit with gaston
	colnames(unadj)[colnames(unadj) == 'pos_Bp'] <- 'pos'
	colnames(unadj)[colnames(unadj) == 'p_value'] <- 'p'
	
	colnames(adj)[colnames(adj) == 'pos_Bp'] <- 'pos'
	colnames(adj)[colnames(adj) == 'p_value'] <- 'p'
	
	treshold = -log10(0.05/dim(x@submaps_list[[1]])[2])  #treshold based on the number of markers
	
	if (save == FALSE) { # Default, just print plots on different windows

	opar <- par(mfrow = c(2,2))  
	
	#dev.new(height = 5.1, width = 5)
	print(p1)
	
	#dev.new(height = 5.1, width = 5, ypos = 650)
	print(p2)
	
	#dev.new(width = 14.25, height= 5.1, ypos = 0, xpos = 640)
	gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"))
	abline(h=treshold, col = 'red')
	
	# adjusted
	#dev.new(width = 14.25, height= 5.1, ypos = 650, xpos = 640)
	gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"))
	abline(h=treshold, col = 'red')
	
	par(opar)
	
	} else { # save plots in png files 
	
	png(paste('qqplot.GLM.', expl_var,'unadj.png')) 		# U = unadjusted
	print(p1)
	dev.off()
	
	png(paste('qqplot.GLM.', expl_var, 'adj.png'))		# A = adjusted
	print(p2)
	dev.off()
	
	png(paste( 'manhattanplot.GLM.', expl_var,'unadj.png'))
	gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"))
	abline(h=treshold, col = 'red')
	dev.off()	
	
	png(paste( 'manhattanplot.GLM.', expl_var,'adj.png'))
	gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"))
	abline(h=treshold, col = 'red')
	dev.off()	
	
	}
			  
		
}



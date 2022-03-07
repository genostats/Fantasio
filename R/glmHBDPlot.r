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
	
	if (save == FALSE) { # Default, just print plots on different windows

	dev.new(height = 5.1, width = 5)
	print(p1)
	
	dev.new(height = 5.1, width = 5, ypos = 650)
	print(p2)
	
	# Manhattan Plot
	# Change colnames to fit with gaston
	colnames(unadj)[colnames(unadj) == 'pos_Bp'] <- 'pos'
	colnames(unadj)[colnames(unadj) == 'p_value'] <- 'p'
	
	dev.new(width = 14.25, height= 5.1, ypos = 0, xpos = 640)
	gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"))
	#abline(h=treshold, col = 'red')
	
	colnames(adj)[colnames(adj) == 'pos_Bp'] <- 'pos'
	colnames(adj)[colnames(adj) == 'p_value'] <- 'p'
	
	# adjusted
	#treshold = -log10(1/dim(x@submaps_list[[1]])[2])  #treshold based on the number of markers
	dev.new(width = 14.25, height= 5.1, ypos = 650, xpos = 640)
	gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"))
	#abline(h=treshold, col = 'red')
	
	
	} else { # save plots in png files 
	
	png(paste('qqplot.GLM.', expl_var,'U.png')) 		# U = unadjusted
	p1
	dev.off()
	
	png(paste('qqplot.GLM.', expl_var, 'A.png'))		# A = adjusted
	p2
	dev.off()
	
	png(paste( 'manhattanplot.GLM.', expl_var,'unadj.png'))
	gaston::manhattan(unadj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- unadjusted data"), chrom.col = c("darksalmon", "darkturquoise"))
	dev.off()	
	
	png(paste( 'manhattanplot.GLM.', expl_var,'adj.png'))
	gaston::manhattan(adj, main = paste("Manhattan Plot \n GLM with ", expl_var, "- adjusted data"), chrom.col = c("darksalmon", "darkturquoise"))
	dev.off()	
	
	}
			  
		
}



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
#' @seealso glmHBDq()
#' 
#' @export

	## Com Margot : à bouger :
	library(qqman)

glmHBDPlot = function ( x, expl_var, covar_df , covar , n.cores = 1, save = FALSE) {
	
	final_unadj 	<- glmHBD(x = x, expl_var = expl_var , n.cores = n.cores)
	message("-----------> GLM on UNADJUSTED data OK \n")

	final 		<- glmHBD(x = x, expl_var = expl_var , covar_df = covar_df , covar = covar, n.cores = n.cores )
	message("-----------> GLM on ADJUSTED data OK \n")
	
	final <- final [which(final$p_value != 0), ]
	final_unadj <- final_unadj [which(final_unadj$p_value != 0), ]
	
	# QQ Plot
	# unadjusted
	p1 <- ggQQPlot(final_unadj$p_value) + ggtitle(paste("QQ-plot GLM with", expl_var, "- unadjusted data"))
	
	# adjusted
	p2 <- ggQQPlot(final$p_value) + ggtitle(paste("QQ-plot GLM with", expl_var, "- adjusted data"))
	
	if (save == FALSE) { # Default, just print plots on different windows

	dev.new(height = 5.2, width = 5.5, xpos = 470)
	print(p1)
	
	dev.new(height = 5.2, width = 5.5, xpos = 1000)
	print(p2)
	
	# Manhattan Plot
	# adjusted
	dev.new(width = 20, height= 5, ypos = 650)
	man <- manhattan(final, main = paste("Manhattan Plot GLM with ", expl_var) ,
			  chr = "chr",
			  bp = "pos_Bp",
			  snp = "snps",
			  p = "p_value",
			  col = c("darksalmon", "darkturquoise") )
	
	} else { # save plots in png files 
	
	png(paste('qqplot.GLM.', expl_var,'U.png')) 		# U = unadjusted
	p1
	dev.off()
	
	png(paste('qqplot.GLM.', expl_var, 'A.png'))		# A = adjusted
	p2
	dev.off()
	
	png(paste( 'manhattanplot.GLM.', expl_var,'.png'))
	man <- manhattan(final, main = paste("Manhattan Plot GLM with ", expl_var) ,
			  chr = "chr",
			  bp = "pos_Bp",
			  snp = "snps",
			  p = "p_value",
			  col = c("darksalmon", "darkturquoise") )
	dev.off()	
	
	}
			  
		
}



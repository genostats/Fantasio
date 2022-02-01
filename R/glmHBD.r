glmHBD <- function( x, expl_var, covar_df, covar, n.cores = 1) {

	if (expl_var == 'HBD_prob') {
		# Recovery HBD_prob
		hbd <- x@HBD_recap
	} else if (expl_var == 'FLOD') {
		# Recovery FLOD
		hbd <- x@FLOD_recap
	} else {
		stop("Explanatory variable must be 'HBD_prob' or 'FLOD'")
	}
		
	# Recovery phenotype
	id <- sub("^\\d*:", "" , row.names(hbd))
	id.index <- match ( id, x@bedmatrix@ped$id )
	pheno <- x@bedmatrix@ped$pheno [id.index]
	pheno <- ifelse( pheno == 2, 1, 0) # Translate phenotype
	
	# Recovery chr, snps, pos_cM and pos_Bp 
	final <- getPosition(x, hbd)	
	
	if(n.cores == 1 ) {
		# unadjusted 
		if (missing(covar_df)) {
			for(i in 1:ncol(hbd)){
				model <- glm( pheno ~ hbd[,i] , family = binomial)
				final[i,'estimate'] 	<- summary(model)$coef[2,1]
				final[i,'std_error'] 	<- summary(model)$coef[2,2]
				final[i,'z_value'] 	<- summary(model)$coef[2,3]
				final[i,'p_value'] 	<- summary(model)$coef[2,4]
			}
		}
		
		# adjusted 
		else {	
			if(missing(covar)) {
				df <- covar_df				 # take all covar given in the dataframe
			} else {
				df <- as.data.frame(covar_df[ , covar]) #rownames covar_df  = individual id 	
			}
			
			for(i in 1:ncol(hbd)){
				model <- glm( pheno ~ hbd[,i] + . , data = df,  family = binomial)
				final[i,'estimate'] 	<- summary(model)$coef[2,1]
				final[i,'std_error'] 	<- summary(model)$coef[2,2]
				final[i,'z_value'] 	<- summary(model)$coef[2,3]
				final[i,'p_value'] 	<- summary(model)$coef[2,4]
			}
		}
	}
	
	# Parallelization - n.cores > 1
	else {
 		## Com Margot : à bouger :
 		library(foreach)
 		cl <- parallel::makeCluster(n.cores)
		doParallel::registerDoParallel(cl)
		
		if(missing(covar_df)) {

			res <- foreach (i = 1:ncol(hbd), .combine = rbind) %dopar% {
				model <- glm( pheno ~ hbd[,i] ,  family = binomial)
				estimate <- summary(model)$coef[2,1]
				std_error <- summary(model)$coef[2,2]
				z_value <- summary(model)$coef[2,3]
				p_value <- summary(model)$coef[2,4]
				data.frame( estimate, std_error, z_value, p_value)
	 		}
	 	}
	 	
	 	else {
	 	
			if(missing(covar)) {
				df <- covar_df				 # take all covar given in the dataframe
			} else {
				df <- as.data.frame(covar_df[ , covar]) #rownames covar_df  = individual id 	
			} 
	 		
			res <- foreach (i = 1:ncol(hbd), .combine = rbind) %dopar% {
				model <- glm( pheno ~ hbd[,i] + . , data = df, family = binomial)
				estimate <- summary(model)$coef[2,1]
				std_error <- summary(model)$coef[2,2]
				z_value <- summary(model)$coef[2,3]
				p_value <- summary(model)$coef[2,4]
				data.frame( estimate, std_error, z_value, p_value)
	 		}
	 	}
	 	
	 	parallel::stopCluster(cl)
	 	final <- cbind(final, res)
	}
	return (final)
} 
	


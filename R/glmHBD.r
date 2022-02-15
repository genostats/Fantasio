#' Logistic regression on HBD prob or FLOD
#' 
#' @param x the bedmatrix after Fantasio()
#' @param expl_var the explanatory variable 'FLOD' or 'HBD_prob'
#' @param covar_df a dataframe containing covariates
#' @param covar covariates of interest 
#' if missing, all covariates of the dataframe are considered
#' @param n.cores number of cores for parallelization calculation (default = 1)
#' @param save choose if plot are saved or not (default = FALSE)
#' @param run whether the fonction is called or not
#' 
#' @export

glmHBD <- function( x, expl_var, covar_df, covar, n.cores = 1, run, phen.code) {
  
  if(class(x)[1] != "atlas")
    stop("Need an atlas")

  if (run) { 
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
	  if (phen.code == 'plink') {
	    pheno <- ifelse(pheno == 1, 0, ifelse(pheno == 2, 1, NA))# Translate phenotype
	  }
	
	  # Recovery chr, snps, pos_cM and pos_Bp 
	  final <- getPosition(x, hbd)	
	
	  if(n.cores == 1 ) {
		  # unadjusted 
		  if (missing(covar_df)) {
		    message("No covariates given for the analysis = unadjusted data. To use covariates import a dataframe.")
			  for(i in 1:ncol(hbd)){
				  model <- glm( pheno ~ hbd[,i] , family = binomial)
				  final[i,'estimate'] 	<- summary(model)$coef[2,1]
				  final[i,'std_error'] 	<- summary(model)$coef[2,2]
				  final[i,'z_value'] 	<- summary(model)$coef[2,3]
				  final[i,'p_value'] 	<- summary(model)$coef[2,4]
			  }
		  x@logisticRegression$unadj <- cbind(final, res)
		  }
		
		  # adjusted 
		  else {	
			  if(missing(covar)) {
			    message("No covariates specified - All covariates of the dataframe will be used.")
				  df <- covar_df				 # take all covar given in the dataframe
			  } else {
				  df <- as.data.frame(covar_df[ id.index , covar]) #rownames covar_df  = individual id 	
			  }
			
			  for(i in 1:ncol(hbd)){
				  model <- glm( pheno ~ hbd[,i] + . , data = df,  family = binomial)
				  final[i,'estimate'] 	<- summary(model)$coef[2,1]
				  final[i,'std_error'] 	<- summary(model)$coef[2,2]
				  final[i,'z_value'] 	<- summary(model)$coef[2,3]
				  final[i,'p_value'] 	<- summary(model)$coef[2,4]
			  }
		  x@logisticRegression$adj <- cbind(final, res)  
		  }
	  }
	
	  # Parallelization - n.cores > 1
	  else {
 		  ## Com Margot : à bouger :
 		  library(foreach)
 		  cl <- parallel::makeCluster(n.cores)
		  doParallel::registerDoParallel(cl)
		
		  if(missing(covar_df)) {
		    message("No covariates given for the analysis = unadjusted data. To use covariates import a dataframe.")
		  
			  res <- foreach (i = 1:ncol(hbd), .combine = rbind) %dopar% {
				  model <- glm( pheno ~ hbd[,i] ,  family = binomial)
				  estimate <- summary(model)$coef[2,1]
				  std_error <- summary(model)$coef[2,2]
				  z_value <- summary(model)$coef[2,3]
				  p_value <- summary(model)$coef[2,4]
				  data.frame( estimate, std_error, z_value, p_value)
			  }
			x@logisticRegression$unadj <- cbind(final, res)  
      }
	 	
	 	  else {
	 	
			  if(missing(covar)) {
			    message("No covariates specified - All covariates of the dataframe will be used.")
				  df <- covar_df				 # take all covar given in the dataframe
			  } else {
				  df <- as.data.frame(covar_df[ id.index , covar]) #rownames covar_df  = individual id 	
			  } 
	 		
			  res <- foreach (i = 1:ncol(hbd), .combine = rbind) %dopar% {
				  model <- glm( pheno ~ hbd[,i] + . , data = df, family = binomial)
				  estimate <- summary(model)$coef[2,1]
				  std_error <- summary(model)$coef[2,2]
				  z_value <- summary(model)$coef[2,3]
				  p_value <- summary(model)$coef[2,4]
				  data.frame( estimate, std_error, z_value, p_value)
			  }
			x@logisticRegression$adj <- cbind(final, res)  
	 	  }
	 	 parallel::stopCluster(cl)
	  }
  }
  x
} 
	


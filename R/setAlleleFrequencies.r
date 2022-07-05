#' This function assign allele frequencies from a reference data set
#' Use this function if you want to change the allelic frequencies in the bed.matrix object.

#' @param x the bed.matrix to study
#' @param x_ref the reference bed.matrix to set allelic frequencies

#' @export

setAlleleFrequencies <- function (x, x_ref)
{
      # look if missing SNPs in the reference dataset
	res <- match(x@snps$id, x_ref@snps$id)
	snps_eliminated <- is.na(res)
	
	if (sum(snps_eliminated) > 0) {
	warning( sum(snps_eliminated), " SNPs are eliminated from the analysis" )
	}
	
      # retrieve indexes of SNPs that match
	x_ind <- which(!is.na(res))
	x_ref_ind <- res[!is.na(res)]
	
      # force the same position for the SNP match
	x_ref@snps$pos[x_ref_ind] <- x@snps$pos[x_ind]				
#### MD_com : ici on travaille avec des données de génotypage. Attention pour les données de séquençage.
	m <- SNP.match(x, x_ref, by = "chr:pos:alleles")			
	
## on ne peut pas éliminer m$index = NA
#### MD_com : (prévoir dans le code C++ que si p = NA on ne tient pas compte du SNP)
	
      # assign allelic frequencies from the reference dataset :
      ## Set NA.values to SNPs eliminated from the study
      ## change value if "swap" (= a SNP A/C would match a SNP C/A)
	x@p <- ifelse(is.na(res), NA, ifelse(m$swap, x_ref@p[match(x@snps$id, x_ref@snps$id)], x_ref@p[match(x@snps$id, x_ref@snps$id)]))
		
      ## avoid extreme values ​​0 and 1 for probability calculations
	x@p[x_ind] <- ifelse(x_ref@p[x_ref_ind] == 1, 1 - 1e-3, x_ref@p[x_ref_ind])
	x@p[x_ind] <- ifelse(x_ref@p[x_ref_ind] == 0, 1e-3, x_ref@p[x_ref_ind])
	
	x
	#rm(x_ind et x_ref_ind) pour libérer de l'espace ? 
}


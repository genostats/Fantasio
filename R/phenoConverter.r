#' Fonction to convert phenotype when the phenotype coded with plink
#' convert plink to R coding phenotype

#' @param x the bedmatrix
#' @param phen.code type of coding phenotype 'R' or 'plink'
#' R : 0:control ; 1:cases ; NA:unknown (default)
#' plink : 1:control ; 2:cases ; 0/-9/NA:unknown

phenoConverter = function (x, phen.code = 'R') {
	if (phen.code == 'plink') {
		x@ped$pheno <- ifelse(x@ped$pheno == 1, 0, ifelse(x@ped$pheno == 2, 1, NA))
	} else continue
	}
	

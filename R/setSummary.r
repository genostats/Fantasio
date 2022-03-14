#' Creation summary statistic files
#' 
#' This function is uses to ouput all the summary files after the creation of the atlas by the following two functions : 
#' `makeAtlasByHotspots` and `makeAtlasByDistance` is over.
#' 
#' 
#' @param atlas a atlas object
#' @param list.id you can either :
#'     - ignore this parameter if you want to compute HBD, FLOD and HFLOD 
#'       for individuals who are considerated inbred and with a quality
#'       greater or equal to 95%}
#'     - enter a list of individual for a computation of HBD, FLOD score HFLOD score for them
#'     - use "all" for a computation of HBD, FLOD score and HFLOD score for every individual
#' @param probs a flag indicating if the HBD probabilities and FLOD score has been computed for the submaps (default is TRUE)
#' @param recap.by.segments if the summary files has to be computed considering segments or snps (defaut is FALSE) 
#' (more information in the documentation of `makeAtlasByHotspots` or`makeAtlasByDistance` functions)
#' @param HBD.threshold value of the HBD probability threshold used to determine whether a segment is HBD or not (default is 0.5)
#' @param q  Allows the user to choose the assumed frequency of the mutation involved in the disease for each individual (default is 0.0001)
#' @param quality The minimum percentage use to assume if a submap is valid (default is 95)
#' @param n.consecutive.markers the number of consecutive markers with a probability equal or greater to the value of the threshold, to be used to find HBDsegments (default is 5)
#' 
#' @details the function add content to the atlas object.
#' @details - summary for likelihood0/likelihood1
#' @details - summary for estimation of a and f
#' @details - summary for the submaps
#' @details - summary for HBD probabilities
#' @details - summary for HFLOD
#' @details - comutation of HBDsegments
#' 
#' @return return an atlas containing every summary created
#' 
#' @seealso Fantasio
#' @seealso makeAtlasByDistance
#' @seealso segmentsListByHotspots
#' @seealso festim
#' @seealso setHBDProb
#' @seealso setFLOD
#' @seealso setHFLOD
#' @seealso recap
#' @seealso submapLikelihood
#' @seealso submapEstim
#' @seealso submapSummary
#' 
#' 
#' @export
setSummary <- function (atlas, list.id, probs = TRUE, recap.by.segments = FALSE,
    q = 1e-04, HBD.threshold = 0.5, quality = 95, n.consecutive.markers = 5, phen.code = c('plink', 'R')) {

  phen.code <- match.arg(phen.code)
  
  if(class(atlas)[1] != "atlas")
    stop("Need an atlas")
 
  atlas@likelihood_summary <- submapLikelihood(atlas@submaps_list)
  atlas@estimation_summary <- submapEstim(atlas@submaps_list)
  atlas@submap_summary <- suppressWarnings(submapSummary(atlas@submaps_list))
  
  if(probs) {
    if (phen.code == 'plink') {
      test <- any( atlas@submap_summary$pheno == 2 ) 
    } else {
      test <- any( atlas@submap_summary$pheno == 1 ) 
    }
    if(missing(list.id)) { # pas de list.id : défaut 
      if(test) { # il y a des atteints
        # on calcule les probas HBD et les FLOD sur les individus consanguins avec qualité suffisante
        # le HFLOD sur les atteints parmi ceux là
        w.HBD   <- which( atlas@submap_summary$quality >= quality & atlas@submap_summary$inbred )
        if (phen.code == 'plink') {
          w.HFLOD <- match( which(atlas@submap_summary$quality >= quality & atlas@submap_summary$inbred & atlas@submap_summary$pheno == 2), w.HBD )
        } else {
          w.HFLOD <- match( which(atlas@submap_summary$quality >= quality & atlas@submap_summary$inbred & atlas@submap_summary$pheno == 1), w.HBD )
        }
      } else {
        # on calcule les probas HBD, les FLOD et les HFLOD sur tous les consanguins 
        # avec qualité
        w.HBD   <- which( atlas@submap_summary$quality >= quality & atlas@submap_summary$inbred )
        w.HFLOD <- seq_along(w.HBD)
        if (phen.code == 'plink') {
          warning("No individual with pheno = 2.\nUsing all inbred individuals with good estimation quality.")
        } else {
          warning("No individual with pheno = 1.\nUsing all inbred individuals with good estimation quality.")
        }
      }
    } else { # on calcule sur les individus donnés !
      if(list.id == "all") {
        w.HBD <- seq_len(nrow(atlas@submap_summary))
        w.HFLOD <- seq_along(w.HBD)
      } else {
        w.HBD <- match( list.id, uniqueIds(atlas@submap_summary$famid, atlas@submap_summary$id) )
        w.HFLOD <- seq_along(w.HBD)
      }
    }

    atlas <- setHBDProbAndFLOD(atlas, w.id = w.HBD, q=q)
    #atlas <- setHBDProb(atlas, w.id = w.HBD)
    #atlas <- setFLOD(atlas, w.HBD, q=q)
   
    # test class of first submap to check if recap is ok 
    if(class(atlas@submaps_list[[1]])[1] == "snpsMatrix" & recap.by.segments) {
      warning("Submaps created by Distance force 'recap.by.segments = FALSE'")
      recap.by.segments <- FALSE
    }
    atlas@bySegments <- recap.by.segments

    # recapitulation !
    #l2 <- recap(atlas, recap.by.segments)
    #atlas@HBD_recap <- l2[[1]]
    #atlas@FLOD_recap <- l2[[2]]

    #atlas@HBDsegments <- HBDSegments(atlas, threshold = HBD.threshold, n.consecutive.markers = n.consecutive.markers) 
    
    #atlas <- setHFLOD(atlas, w.HFLOD)
  }
  atlas
}


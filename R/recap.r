##################################################################################
#This function creates HBD and FLOD recap dataframe                              #
#                                                                                #
#!!! submaps : the list of object                                                #                                       
#!!! by_segments : whether you want the recap by segments or not                 #
#!!! list.id : a list individual                                                 #
#                                                                                # 
#*** return a list of two dataframes with HBD and FLOD recap in it               #
##################################################################################

recap <- function(submaps, by_segments=F, list.id)
{
  
  if(!missing(list.id))
  {
    if(list.id == "all")
    {
      condition <- 1:nrow(submaps@submap_summary)
    }else{
      vec <- strsplit(list.id, "_")
      condition <- sapply(vec, function(i) which(submaps@submap_summary$FID == i[1] & submaps@submap_summary$IID == i[2]))
    }
  }else{
    condition <- which(submaps@submap_summary$QUALITY >= 95 & submaps@submap_summary$INBRED)
  }
  
  if(length(submaps@atlas) == 1)
  {
  	matrice_HBD <- submaps@atlas[[1]]@HBD.prob
  	matrice_FLOD <- submaps@atlas[[1]]@FLOD
  	l <- list(matrice_HBD, matrice_FLOD)
  	return(l)	
  }
  
  #recuperer les probas HBD pour chaque sous cartes sous forme d'une liste
  proba_HBD  <- submap.HBD(submaps)
  proba_FLOD <- submap.FLOD(submaps)
  
  if(by_segments)
  {
    recap.by.segments(submaps, proba_HBD, proba_FLOD)
  }
  else
  {
    recap.by.snps(submaps, proba_HBD, proba_FLOD)
  }
}  

##################################################################################
#This function creates a dataframe with HBD.segments used for plottting graphs   #
#by segments                                                                     #
#                                                                                #
#!!! Submaps : the list of objects                                               #                                       
#!!! HBD_recap : a dataframe with HBD probabilities                              #
#!!! n.consecutive.marker : the number of consecutive HBD probabilities >        #
#    threshold                                                                   #
#                                                                                #
#*** return a dataframe with HBD segments                                        #
##################################################################################

HBD.segments.by.segments <- function(submaps, HBD_recap, n.consecutive.marker, threshold)
{
  l <- list()#liste des individus
  
  individuals_name <- rownames(HBD_recap)#get the name of the individual
  individuals_name <- strsplit(individuals_name, "_")
  individuals_name <- sapply(individuals_name, function(i) match(i, submaps@bedmatrix@ped$id))
  individuals_name <- individuals_name[!is.na(individuals_name)]
  
  #find the status of the individual
  STATUS <- submaps@bedmatrix@ped$pheno[individuals_name]
  family_id <- submaps@bedmatrix@ped$famid[individuals_name]
  individuals_name <- submaps@bedmatrix@ped$id[individuals_name]
  
  #recuperer le nom des marqueurs pour connaitre l'emplacement du segments dans le chr
  marker_names <- colnames(submaps@atlas[[1]]@HBD.prob)
  
  
  #a quel chromosome appartient ce marqueur
  correspondance <- match(marker_names, submaps@bedmatrix@snps$id)
  
  chromosome <- submaps@bedmatrix@snps$chr[correspondance]
  
  ##recuperer le nom de l'individus
  #nom <- rownames(HBD_recap)
  #nom <- strsplit(nom, "_")
  
  #recuperer les distances et positions de chaques segment 
  
  start_dist <- c()
  end_dist <- c()
  
  start_pos <- c() 
  end_pos <- c()
  
  for(i in 1:22)#boucle sur les chr
  {
    segment_chr <- submaps@segments_list[[i]]#recuperer la liste des segments pour le chr i
    
    for(j in 1:length(segment_chr))#boucle sur les segments du chromosome
    {
      if(length(segment_chr[[j]])== 0) next()#si pas de segment donc on passe
      
      start <- segment_chr[[ j ]][ 1 ]# premier marqueur
      end <- segment_chr[[ j ]][ length( segment_chr[[ j ]] ) ]# dernier marqueur
      
      #debut du segment 
      start_dist <- c(start_dist, submaps@bedmatrix@snps$dist[ start ])
      
      #fin du segment -> dernier marqueur
      end_dist   <- c(end_dist, submaps@bedmatrix@snps$dist[ end ] )
      
      #debut du segment 
      start_pos  <- c(start_pos, submaps@bedmatrix@snps$pos[ start ])
      
      #fin du segment -> dernier marqueur 
      end_pos    <- c(end_pos, submaps@bedmatrix@snps$pos[ end ])
      
      
    }
  }
  
  
  
  #create the dataframe for every individual
  for(i in 1:nrow(HBD_recap))
  {
    data<-c(HBD_recap[i,])#ligne i de la matrice HBD_recap
    test<- as.vector(data >= threshold)#test
    
    min_segment_size = n.consecutive.marker
    
    all_segments<-rle(test)
    
    good_segments<-as.numeric( which( all_segments$length >= min_segment_size & all_segments$value ) )
    
    good_segments_length <-as.numeric( all_segments$length[ good_segments ] )
    
    #good_segments_start<- as.numeric( ( c(0,cumsum(all_segments$lengths))+1)[ good_segments ] )
    
    good_segments_start<- as.numeric( cumsum(all_segments$lengths)[ good_segments ] )
    
    good_segments_end<-as.numeric( good_segments_start+good_segments_length )
    
    
    
    #splitting <- strsplit(nom[i], "_")#recuperer l'id de l'individus
    
    
    segment_dataframe<-data.frame(individual  = rep(individuals_name[i], length(good_segments_start)),
                                  family      = rep(family_id[i], length(good_segments_start)),
                                  status      = rep(STATUS[i], length(good_segments_start)),
                                  start=good_segments_start, 
                                  end=good_segments_end,
                                  size=good_segments_length,
                                  chromosome  = chromosome[good_segments_start],
                                  start_pos   = start_pos[good_segments_start],
                                  end_pos     = end_pos[good_segments_end],
                                  start_dist  = start_dist[good_segments_start],
                                  end_dist    = end_dist[good_segments_end])
    
    #treating the case when segments overlaps two different chromosomes
    overlap <- which(segment_dataframe$start_dist > segment_dataframe$end_dist)
    
    
    #treating the case when segments overlaps two different chromosomes
    if(length(overlap) != 0)
    {
      segment_dataframe <- segment_dataframe[-overlap,,drop=F]
    }
    
    l[[i]] <- segment_dataframe
  }
  
  return(l)
}


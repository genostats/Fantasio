HBD.segments <- function(submaps, HBD.recap, n.consecutive.marker = 5, threshold = 0.5, by_segments=F)
{
  if(by_segments)
  {
    l <- list()
    marker_names <- colnames(submaps@atlas[[1]]@HBD.prob)#recuperer le nom des marqueurs pour connaitre l'emplacement du segments dans le chr
    chromosome <- submaps@bedmatrix@snps$chr[match(marker_names, submaps@bedmatrix@snps$id)]#a quel chromosome appartient ce marqueur
    
    nom <- rownames(HBD.recap)#get the name of the individual
    
    
    #get distance and position for every segments
    
    start_dist <- c()
    end_dist <- c()
    start_pos <- c() 
    end_pos <- c()
    
    for(i in 1:22)#boucle sur les chr
    {
      segment_chr <- submaps@segments_list[[i]]#recuperer la liste des segments pour le chr i
      for(j in 1:length(segment_chr))#boucle sur les segments du chromosome
      {
        #if(length(x@snps$dist[ segment_chr[[ j ]][ 1 ] ] == 0)) next()#si une valeur vide pas de segment donc on passe
        
        start_dist <- c(start_dist, submaps@bedmatrix@snps$dist[ segment_chr[[ j ]][ 1 ] ])#debut du segment
        end_dist   <- c(end_dist, submaps@bedmatrix@snps$dist[ segment_chr[[ j ]][ length( segment_chr[[ j ]] ) ] ] )#fin du segment
          
        start_pos  <- c(start_pos, submaps@bedmatrix@snps$pos[ segment_chr[[ j ]][ 1 ] ])#debut du segment
        end_pos    <- c(end_pos, submaps@bedmatrix@snps$pos[ segment_chr[[ j ]][ length( segment_chr[[ j ]] )] ])#fin du segment
        
        
      }
    }
    
    #create the dataframe for every individual
    for(i in 1:nrow(HBD.recap))
    {
      data<-c(HBD.recap[i,])#save the line 
      test<- (data >= threshold)#test
      
      min_segment_size = n.consecutive.marker
      all_segments<-rle(test)
      
      good_segments<-which( all_segments$length >= min_segment_size & all_segments$value )
      
      good_segments_length <-all_segments$length[ good_segments ]
      
      good_segments_start<- ( c(0,cumsum(all_segments$lengths)) +1)[ good_segments ]
      
      good_segments_end<-good_segments_start+good_segments_length-1
      
      splitting <- strsplit(nom[i], "_")#recuperer l'id de l'individus
      
      STATUS <- submaps@bedmatrix@ped$pheno[which(submaps@bedmatrix@ped$id == splitting[[1]][1])]#recuperer son status
  
      segment_dataframe<-data.frame(individual  = rep(nom[i], length(good_segments)),
                                    status      = rep(STATUS, length(good_segments)),
                                    start=good_segments_start, 
                                    end=good_segments_end,
                                    size=good_segments_length,
                                    chromosome  = chromosome[as.numeric(good_segments_start)],
                                    start_pos   = start_pos[as.numeric(good_segments_start)],
                                    end_pos     = end_pos[as.numeric(good_segments_end)],
                                    start_dist  = start_dist[good_segments_start],
                                    end_dist    = end_dist[as.numeric(good_segments_end)])
      l[[i]] <- segment_dataframe
    }
    return(l)
  }
  
  
#######################################################################################
#######################################################################################
######################################################################################
  
  else
  {
    l <- list()
    nom <- rownames(HBD.recap)#get the name of the individual
    
    for(i in 1:nrow(HBD.recap))
    {
      data<-c(HBD.recap[i,])#save the line 
      marker <- colnames(HBD.recap)#save the marker
      test<- (data >= threshold)#test
      
      correspondance <- match(colnames(HBD.recap), submaps@bedmatrix@snps$id)#match between marker's name in HBD.recap and the bedmatrix
      chr <- submaps@bedmatrix@snps$chr[correspondance]#chromosome on which we have the marker
      
      min_segment_size = n.consecutive.marker#minimum size of the marker
      
      # get the segments
      
      all_segments<-rle(test)
      
      good_segments<-which( all_segments$length >= min_segment_size & all_segments$value )
      
      good_segments_length <-all_segments$length[ good_segments ]
      
      good_segments_start<- (c(0,cumsum(all_segments$lengths)) +1)[ good_segments-1 ]
      
      good_segments_end<-good_segments_start+good_segments_length-1
      
      #chromosome delimitation
      chr_segments_start <- chr[as.numeric(good_segments_start)]
      
      chr_segments_end <- chr[as.numeric(good_segments_end)]
      
      overlap <- which(chr_segments_start != chr_segments_end)
      
      #treating the case when segments overlaps two different chromosomes
      if(length(overlap) != 0)
      {
        #put to NA all the values
        chr[overlap] <- NA
        good_segments_start[overlap] <- NA
        good_segments_end[overlap] <- NA
        good_segments_length[overlap] <- NA
        chr_segments_start[overlap] <- NA
        
        #removing NA value
        
        chr <- chr[!is.na(chr)]
        good_segments_start <- good_segments_start[!is.na(good_segments_start)]
        good_segments_end <- good_segments_end[!is.na(good_segments_end)]
        good_segments_length <- good_segments_length[!is.na(good_segments_length)]
        chr_segments_start <- chr_segments_start[!is.na(chr_segments_start)]
      }
      
      #finding distance and position 
      
      start_pos <- submaps@bedmatrix@snps$pos[correspondance[as.numeric(good_segments_start)]]
      end_pos <- submaps@bedmatrix@snps$pos[correspondance[as.numeric(good_segments_end)]]
      
      
      start_dist <- submaps@bedmatrix@snps$dist[correspondance[as.numeric(good_segments_start)]]
      end_dist <- submaps@bedmatrix@snps$dist[correspondance[as.numeric(good_segments_end)]]
      
      
      
      
      #find the status of the individual
      splitting <- strsplit(nom[i], "_")
      STATUS <- submaps@bedmatrix@ped$pheno[which(submaps@bedmatrix@ped$id == splitting[[1]][1])]
      
      
      
      #dataframe
      
      segment_dataframe<-data.frame(individual = rep(nom[i], length(start_pos)),
                                    status = rep(STATUS, length(start_pos)),
                                    start=good_segments_start, 
                                    end=good_segments_end ,
                                    size=good_segments_length,
                                    chromosome=chr_segments_start,
                                    start_pos = start_pos,
                                    end_pos = end_pos,
                                    start_dist = start_dist,
                                    end_dist = end_dist)
      
      
      
      l[[i]] <- segment_dataframe
    }
    return(l)
  }
  
  
  

}


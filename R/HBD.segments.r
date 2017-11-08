HBD.segments <- function(submaps, HBD.matrix, n.consecutive.marker = 5, threshold = 0.5)
{
  l <- list()
  nom <- rownames(HBD.matrix)#get the name of the individual
  
  for(i in 1:nrow(HBD.matrix))
  {
    data<-c(HBD.matrix[i,])#save the line 
    marker <- colnames(HBD.matrix)#save the marker
    test<- (data >= threshold)#test
  
    correspondance <- match(colnames(HBD.matrix), x@snps$id)#match between marker's name in HBD.matrix and the bedmatrix
    chr <- x@snps$chr[correspondance]#chromosome on which we have the marker
    
    min_segment_size = n.consecutive.marker#minimum size of the marker
    
    # get the segments
    
    all_segments<-rle(test)
    
    good_segments<-which( all_segments$length >= min_segment_size & all_segments$value )
    
    good_segments_length <-all_segments$length[ good_segments ]
    
    good_segments_start<-( cumsum(all_segments$lengths)+1 )[good_segments-1]
    
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
    
    start_pos <- x@snps$pos[correspondance[as.numeric(good_segments_start)]]
    end_pos <- x@snps$pos[correspondance[as.numeric(good_segments_end)]]
    
    
    start_dist <- x@snps$dist[correspondance[as.numeric(good_segments_start)]]
    end_dist <- x@snps$dist[correspondance[as.numeric(good_segments_end)]]
    
    
    
    
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
  
  l

}


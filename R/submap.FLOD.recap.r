submap.FLOD <- function(h)
{
  if(class(h@atlas[[1]])[1] != "snps.matrix" & class(h@atlas[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat, in the list.submaps object.") 
  
  #for each marker gives the HBD proba for all the individuals
  # lapin <- sapply(seq_along(h), function(i) h[[i]]@HBD.prob, simplify="array")
  lapin <- lapply(h@atlas, function(hh) hh@FLOD)
  return(lapin)
  
}







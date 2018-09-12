submap.FLOD <- function(h)
{
  if(class(h@atlas[[1]])[1] != "snpsMatrix" & class(h@atlas[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps, in the submapsList object.") 
    
  #for each marker gives the HBD proba for all the individuals
  lapin <- lapply(h@atlas, function(hh) hh@FLOD)
  return(lapin)
  
}







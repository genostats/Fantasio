submap.HBD <- function(h)
{
  if(class(h@submaps_list[[1]])[1] != "snpsMatrix" & class(h@submaps_list[[1]])[1] != "HostspotsMatrix")
    stop("need either an hotspots.segments list of submaps or a snpsSegments list of submaps.") 
  
  #for each marker gives the HBD proba for all the individuals
  lapin <- lapply(h@submaps_list, function(hh) hh@HBD.prob)
  return(lapin)

}







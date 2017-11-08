submap.HBD <- function(h)
{
  
  #for each marker gives the HBD proba for all the individuals
  # lapin <- sapply(seq_along(h), function(i) h[[i]]@HBD.prob, simplify="array")
  lapin <- sapply(h@atlas, function(hh) hh@HBD.prob, simplify="array")
  return(lapin)

}







HBD.submap.info <- function(h)
{
  lapin <- sapply(seq_along(h), function(i) h[[i]]@HBD.prob, simplify="array")
  return(lapin)
  #sapply(seq_along(h), function(i) h[[i]]@f )
  #rowMeans( lapin[148, , ])
  #sapply(seq_along(h), function(i) h[[i]]@f )[148,]
}



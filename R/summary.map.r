summary.map <- function(h)
{
  snps <- as.vector(sapply(h@atlas, function(x) colnames(x@HBD.prob))) #every markers that has been choosen on submaps
  b <- as.data.frame(table(snps),  stringsAsFactors=FALSE)
  b
}


summary.map <- function(submaps)
{
  snps <- as.vector(sapply(submaps, function(x) colnames(x@HBD.prob))) #every markers that has been choosen on submaps
  b <- as.data.frame(table(snps),  stringsAsFactors=FALSE)
  b
}


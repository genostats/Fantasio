summaryMap <- function(submaps)
{
  snps <- unlist(sapply(submaps, function(x) x@map$id)) #every markers that has been choosen on submaps
  b <- as.data.frame(table(snps),  stringsAsFactors=FALSE)
  b
}


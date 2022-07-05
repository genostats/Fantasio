nullRemover <- function(list)
{
  remove <- sapply(list, function(l) is.null(l))
  list <- list[!remove]
  list
}

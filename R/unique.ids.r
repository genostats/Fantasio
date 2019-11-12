unique.ids <- function(famid, id) paste(famid, id, sep = ":") 

get.famid <- function(x) {
  sapply( strsplit(x, ":"), function(x) x[1] )
}

get.id <- function(x) {
  sapply( strsplit(x, ":"), function(x) x[2] )
}

cleanHotspots <- function(l, minMarkers) {
  index <- c()
  i <- 1L
  while(i <= length(l)) {
    if(1L + diff(l[[i]]) < minMarkers) {
      if(i == length(l)) { # dernier segment
        # on fusionne (i-1) et i, on vire i, on termine
        l[[i-1L]] <- c(l[[i-1L]][1L], l[[i]][2L])
        l[[i]] <- NULL
        return(l)
      } else {
        # on fusionne i et i+1,  on vire i
        l[[i+1L]] <- c(l[[i]][1L], l[[i+1L]][2L])
        l[[i]] <- NULL
      }
    } else { # on n'incrémente i que si il n'y a pas eu fusion de i et i+1 !
      i <- i + 1L
    }
  }
  l
}


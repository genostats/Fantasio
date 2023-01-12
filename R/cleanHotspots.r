cleanHotspots <- function(l, minMarkers) {
  index <- c()
  i <- 1L
  c <- 0
  while(i <= length(l)) {
    if(length(l[[i]])==1) {
      if(1 < minMarkers) {
        if(i == length(l)) { # dernier segment
          # on fusionne (i-1) et i, on vire i, on termine
          l[[i-1L]] <- c(l[[i-1L]][1L], l[[i]][1L])
          l[[i]] <- NULL
          return(l)
        } else {
          # on fusionne i et i+1,  on vire i
          if(length(l[[i+1]]==1)) {
            l[[i+1L]] <- c(l[[i]][1L], l[[i+1L]][1L])
            l[[i]] <- NULL
          } else {
            l[[i+1L]] <- c(l[[i]][1L], l[[i+1L]][2L])
            l[[i]] <- NULL            
          }
        }
      } else {
        i <- i + 1L
      }
    } else {
      if(1L + diff(l[[i]]) < minMarkers) {
        if(i == length(l)) { # dernier segment
          # on fusionne (i-1) et i, on vire i, on termine
          l[[i-1L]] <- c(l[[i-1L]][1L], l[[i]][2L])
          l[[i]] <- NULL
          return(l)
        } else {
          # on fusionne i et i+1,  on vire i
          if(length(l[[i+1]])==1){
            l[[i+1L]] <- c(l[[i]][1L], l[[i+1L]][1L])
            l[[i]] <- NULL
          } else {
            l[[i+1L]] <- c(l[[i]][1L], l[[i+1L]][2L])
            l[[i]] <- NULL            
          }
        }
      } else {
        i <- i + 1L
      }
    }
  }
l
}

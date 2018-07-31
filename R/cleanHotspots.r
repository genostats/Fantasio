cleanHotspots <- function(l, number_of_marker)
{
  index <- c()
  for(i in 1:length(l))
  {
    if(length(l[[i]]) < number_of_marker)
    {
      if(i == length(l))
        l[[i-1]] <- c(l[[i]], l[[i-1]])
      else
        l[[i+1]] <- c(l[[i]], l[[i+1]])
      }else{
        index <- c(index, i)
      }
  }
  l <- l[index]
  l
}
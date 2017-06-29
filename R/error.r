submap.error <- function(f)
{
  out <- tryCatch(
    {
      readLines(con = f, warn = FALSE)
    },
    
    error <- function(cond)
      {
        message(paste("I need an bedmatrix to work, please feed me with one !
                      You might need to use read.bed.matrix(), 
                      Please go check ?read.bed.matrix help"))
        message("Here's the original error message:")
        message(cond)
        return(NA)
      },
    )
  return(out)
}


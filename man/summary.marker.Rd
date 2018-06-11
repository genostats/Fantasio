\name{summary.marker }
\alias{summary.marker }
\title{Give the number of times and the number of marker picked during the creation of the submaps.}
\description{
  Use this function if you want to get the number of time and number of marker picked.
  }
  \usage{marker <- summary.marker(h, x)}
  \arguments{
    \item{h}{a submap.matrix}
    \item{x}{a bed.matrix}
  }
  
  \value{Return a dataframe with two columns 
    \itemize{
      \item{number_of_time_picked}{the number of time markers has been picked}
      \item{number_of_markers}{the number of marker picked}
    }
  }
  \examples{
    require(FEstim)
    system.file("geb")    #a changer
    x <- read.bed.matrix(geb)
    s <- segments(x)
    h <- submaps(x, 1, s)
    marker <- summary.marker(h, x)
    str(estim)
  }
  
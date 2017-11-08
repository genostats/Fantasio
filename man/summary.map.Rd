\name{summary.map }
\alias{summary.map }
\title{The number of time a marker has been chosen on every submap.}
\description{
  Use this function if you want to get a dataframe with the number of time a marker has been chosen on every submap.
  }
  \usage{map <- summary.map(h)}
  \arguments{
    \item{h}{a submap.matrix object}
  }
  
  \value{Return a dataframe with two columns 
  \itemize{
    \item{Maker}{the marker id}
    \item{Freq}{the number of time the marker has been picked}
  }
  }
  \examples{
    require(FEstim)
    system.file("geb")    #a changer
    x <- read.bed.matrix(geb)
    s <- segments(x)
    h <- submaps(x, 1, s)
    map <- summary.map(h)
    str(estim)
  }
  
\name{submap.estim }
\alias{submap.estim }
\title{Give a list of all the a and f computed in the submap.matrix object}
\description{
  Use this function if you want to get all a and f computed by festim from the submap.matrix object
  }
  \usage{estim <- submap.estim(h)}
  \arguments{
    \item{h}{a submap.matrix object}
  }
  
  \value{return a list of two dataframe, the first contains all the f computed and the second all the a computed. 
}
  \examples{
    require(FEstim)
    system.file("geb")    #a changer
    x <- read.bed.matrix(geb)
    s <- segments(x)
    h <- submaps(x, 1, s)
    estim <- submap.estim(h)
    str(estim)
  }
  
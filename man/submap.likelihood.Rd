\name{submap.likelihood}
\alias{submap.likelihood}
\title{Give all the likelihood from an submap.matrix object}
\description{
  Use this function if you want to get a list of two dataframes of all the likelihood (0 and 1) from an submap.matrix object.
  }
  \usage{likelihood <- submap.likelihood(h)}
  \arguments{
    \item{h}{a submap.matrix object}
  }
  
  \value{return a list of two dataframes with likelihood for every submaps in the submap.matrix object.}
  \examples{
    require(FEstim)
    system.file("geb")    #a changer
    x <- read.bed.matrix(geb)
    s <- segments(x)
    h <- submaps(x, 1, s)
    likelihood <- submap.likelihood(h)
    str(likelihood)
  }
  
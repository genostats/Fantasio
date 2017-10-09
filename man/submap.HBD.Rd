\name{submap.HBD}
\alias{submap.HBD}
\title{Create a list with all the HBD slot from a submap.matrix object}
\description{
  Use this function if you want to get a list of all the HBD slot from a submap.matrix object
  }
  \usage{ proba <- submap.HBD(h)}
  \arguments{
    \item{h}{a submap.matrix object}
  }

  \value{return a list with every HBD slot of the submap.matrix object}
  \examples{
    require(FEstim)
    system.file("geb")    #a changer
    x <- read.bed.matrix(geb)
    s <- segments(x)
    h <- submaps(x, 1, s)
    proba <- submap.HBD(h)
    str(proba)
  }
  
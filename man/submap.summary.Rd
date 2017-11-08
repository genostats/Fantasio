\name{submap.summary }
\alias{submap.summary }
\title{Give a summary of all the submap on the submap.matrix object}
\description{
  Use this function if you want to get a summary about all the submaps on the submap.matrix object.
  }
  \usage{ssummary <- submap.summary(h, a.threshold = 1)}
  \arguments{
    \item{h}{a submap.matrix object}
    \item{a.threshold}{a threshold for "a" values, default is 1}
  }
  
  \value{return a dataframe with mutiple columns : 
  \itemize{
    \item{FID}{family identifier}
    \item{IID}{individual identifier}
    \item{STATUS}{status (1 non-affected, 2 affected, 0 unknown)}
    \item{SUBMAPS}{number of submaps used}
    \item{QUALITY}{percentage of valid submaps (i.e. submaps with a < 1)}
    \item{F_MIN}{minimum f on valid submaps}
    \item{F_MAX}{maximum f on valid submaps}
    \item{F_MEAN}{mean f on valid submaps}
    \item{F_MEDIAN}{median f on valid submaps (recommended to estimate f)}
    \item{A_MEDIAN}{median a on valid submaps (recommended to estimate a)}
    \item{pLRT_MEDIAN}{median p-value of LRT tests on valid submaps}
    \item{INBRED}{a flag indicating if the individual is inbred (pLRT_MEDIAN <0.05) or not}
    \item{pLRT_<0.05}{number of valid submaps with a LRT having a p-value below 0.05}
  }}
  \examples{
    require(FEstim)
    system.file("geb")    #a changer
    x <- read.bed.matrix(geb)
    s <- segments(x)
    h <- submaps(x, 1, s)
    ssummary <- submap.sumamry(h)
    str(ssummary)
  }
  
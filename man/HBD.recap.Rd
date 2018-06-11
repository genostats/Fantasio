\name{HBD.recap }
\alias{HBD.recap }
\title{Creation of a dataframe with all the HBD mean probability computed}
\description{
  Use this function if you want to get a datadrame with all HBD mean probability for every marker} 
  \usage{hbdRecap <- HBD.recap(h, by_segments=F)}
  \arguments{
    \item{h}{a submap.matrix}
    \item{by_segments}{whether to use a marker as a whole segment or not, default is false}
  }
  
  \value{Return a dataframe with HBD mean probability computed for every marker or, is by_segment = true, return a dataframe with HBD mean probability computed for every marker picked in a segment. 

  }
  \details{By row you have individuals.
           By column you have markers }
  \examples{
    require(FEstim)
    system.file("geb")    #a changer
    x <- read.bed.matrix(geb)
    s <- segments(x)
    h <- submaps(x, 1, s)
    hbdRecap <- HBD.recap(h)
    str(estim)
}
  
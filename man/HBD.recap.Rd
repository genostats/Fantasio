\name{HBD.recap }
\alias{HBD.recap }
\title{Creation of a dataframe with all the HBD mean probability computed}
\description{
  Use this function if you want to get a datadrame with all HBD mean probability for every marker} 
  \usage{hbdRecap <- HBD.recap(h)}
  \arguments{
    \item{h}{a submap.matrix}
  }
  
  \value{Return a dataframe with HBD mean probability computed for every marker. 

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
  
\name{marker_summary}
\alias{marker_summary}
\title{ computation of a, f, liklelihood0, likelihood1, p.lrt and if probs = TRUE :  HBD, FLOD and HFLOD probabilities }
\description{ Computation of different probabilities on genomic data 
  (a, f, likelihood0, likelihood0, p.lrt, HBD, FLOD, HFLOD) }
\usage{ 
  x <- festim(x, verbose = TRUE, debug = FALSE, probs = TRUE) 
}

\arguments{
  \item{x}{ an msat.matrix or a submap.matrix }
  \item{verbose = TRUE}{if default is not changed, outputs a message : "Estimation f and a for individual # ...", else nothing}
  \item{debug = FALSE}{if default is not changed, outputs nothing, else outputs computation details about a and f}
  \item{probs = TRUE}{Allows computation of HBD, FLOD and HFLOD probabilities}
}

\value{
  An msat.matrix or an submap.matrix with slot a, f, liklelihood0, likelihood1, p.lrt and 
  if probs = TRUE :  HBD, FLOD and HFLOD probabilities computed
}

\seealso{ \code{\link{set.HBD.prob}}
  \code{\link{set.FLOD}}
  \code{\link{set.HFLOD}}
}

\examples{
  require(FEstim)
  system.file("extdata/in_markall.txt", "extdata/in_datall.txt", package="FEstim")
  x <- read.msat.matrix("Example/in_markall.txt", "Example/in_datall.txt")
  x <- festim(x)
  str(x)
}

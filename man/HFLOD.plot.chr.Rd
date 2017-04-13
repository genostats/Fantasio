\name{HFLOD.plot.chr}
\title{graphical output for a given chromosom's positions in function of their HFLOD}
\description{Make a plot of a given chromosom's positions in function of their HFLOD}
\usage{HFLOD.plot.chr(x, chr)}
\arguments{
  \item{x}{a matrix with HFLOD probabilities already computes}
  \item{chr}{a chromosom which HFLOD plot is wanted}
  }
\seealso{ \code{\link{set.HFLOD}},  \code{\link{HFLOD.plot.chrs.png}},  \code{\link{HFLOD.manhattan.plot}}}
\examples{
require(FEstim)
x <- read.msat.matrix("Example/mapfile.in", "Example/datafile.txt")
x <- festim(x)
x <- HBD.prob(x)
x <- set.FLOD(x)
x <- set.HFLOD(x)
HFLOD.plot.chr(x, chr = 7)
}
\name{HFLOD.plot.chrs.png}
\title{PLot of all chromosom's HFLOD in function of their positions}
  \description{output in a folder called HFLOD all the all chromosom's positions in function of their HFLOD}
\usage{HFLOD.plot.chrs.png(x)}
\arguments{
  \item{x}{a matrix with HFLOD probabilities already computes}
}
\seealso{ \code{\link{set.HFLOD}},  \code{\link{HFLOD.plot.chr}},  \code{\link{HFLOD.manhattan.plot}}}
\examples{
  require(FEstim)
  x <- read.msat.matrix("Example/mapfile.in", "Example/datafile.txt")
  x <- festim(x)
  x <- HBD.prob(x)
  x <- set.FLOD(x)
  x <- set.HFLOD(x)
  HFLOD.plot.chrs.png(x)
}
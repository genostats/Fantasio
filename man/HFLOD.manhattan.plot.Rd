\name{HFLOD.manhattan.plot}
\alias{HFLOD.manhattan.plot}
\title{Manhattan plot}
  \description{output a manhattan plot of all chromosom's positions in function of their HFLOD}
\usage{HFLOD.plot.chrs.png(x)}
\arguments{
  \item{x}{a matrix with HFLOD probabilities already computes}
}
\seealso{ \code{\link{set.HFLOD}},  \code{\link{HFLOD.plot.chr}},  \code{\link{HFLOD.plot.chr}}, \code{\link{festim}}}
\examples{
  require(FEstim)
  system.file("extdata/map_snp.in", "extdata/data_snp.in", package="FEstim")
  x <- read.msat.matrix("map_snp.in", "data_snp.in")
  x <- festim(x)
  x <- HBD.prob(x)
  x <- set.FLOD(x)
  x <- set.HFLOD(x)
  HFLOD.manhattan.plot(x)
}
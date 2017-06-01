\name{HFLOD.plot.chr}
\alias{HFLOD.plot.chr}
\title{Graphical output for a given chromosom}
\description{Make a plot of a given chromosom in function of their HFLOD}
\usage{HFLOD.plot.chr(x, chr)}
\arguments{
  \item{x}{a matrix with HFLOD probabilities already computes}
  \item{chr}{a chromosom which HFLOD plot is wanted}
}
\seealso{ \code{\link{set.HFLOD}},  \code{\link{HFLOD.plot.chrs.png}},  \code{\link{HFLOD.manhattan.plot}}, \code{\link{festim}}}
\examples{
require(FEstim)
system.file("extdata/map_snp.in", "extdata/data_snp.in", package="FEstim")
x <- read.msat.matrix("map_snp.in", "data_snp.in")
x <- festim(x)
x <- HBD.prob(x)
x <- set.FLOD(x)
x <- set.HFLOD(x)
HFLOD.plot.chr(x, chr = 7)
}
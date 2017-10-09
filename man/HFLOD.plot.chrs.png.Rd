\name{HFLOD.plot.chrs.png}
\alias{HFLOD.plot.chrs.png}
\title{PLot of all chromosom}
  \description{Output in your current directory, in a folder called "HFLOD" all the  chromosom's positions in function of their HFLOD in ".png" format}
\usage{HFLOD.plot.chrs.png(x)}
\arguments{
  \item{x}{a matrix with HFLOD probabilities already computes (which means it have been passed to festim() with default value )}
}
\seealso{ \code{\link{set.HFLOD}},  \code{\link{HFLOD.plot.chr}},  \code{\link{HFLOD.manhattan.plot}}, \code{\link{festim}}}
\examples{
  require(FEstim)
  system.file("extdata/map_snp.in", "extdata/data_snp.in", package="FEstim")
  x <- read.msat.matrix("map_snp.in", "data_snp.in")
  x <- festim(x)
  HFLOD.plot.chrs.png(x)
}
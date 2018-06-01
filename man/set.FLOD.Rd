\name{set.FLOD}
\alias{set.FLOD}
\title{Computes FLOD probabilities}
\description{Returns the same object with FLOD computes for each marker of each individual }
\usage{FLOD.prob(x, q)}
\arguments{
  \item{x}{an msat.matrix or submap.matrix with HBD probabilities already computes (which means it have been passed to festim() with default values.}
  \item{q}{the assumed frequency of the mutation involved in the 
  disease for each individual. Default is 0.00001.}
}
\seealso{ \code{\link{set.HFLOD}}, \code{\link{festim}}}
\value{a matrix with FLOD values}
\examples{
require(FEstim)
system.file("extdata/map_snp.in", "extdata/data_snp.in", package="FEstim")
x <- read.msat.matrix("map_snp.in", "data_snp.in")
x <- festim(x, probs = FALSE)
x <- HBD.prob(x)
x <- set.FLOD(x, q = 0.00001)
}
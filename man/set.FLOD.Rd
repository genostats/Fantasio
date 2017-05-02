\name{set.FLOD}
\alias{set.FLOD}
\title{Computes FLOD probabilities}
\description{Returns a matrix with FLOD computes for each marker of each individual }
\usage{FLOD.prob(x, q)}
\arguments{
  \item{x}{a matrix with HBD probabilities already computes (which means it have been passed to festim() and HBD.prob() )}
  \item{q}{the assumed frequency of the mutation involved in the 
  disease for each individual. Default is 0.00001.}
}
\seealso{ \code{\link{set.HFLOD}}, \code{\link{festim}}}
\value{a matrix with FLOD values}
\examples{
require(FEstim)
x <- read.msat.matrix("Example/mapfile.in", "Example/datafile.txt")
x <- festim(x)
x <- HBD.prob(x)
x <- set.FLOD(x, q = 0.00001)
}
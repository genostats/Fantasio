\name{set.FLOD}
\title{Computes FLOD probabilities for all the individual's markers}
\description{Returns a matrix with FLOD computes for each marker of each individual}
\usage{FLOD.prob(x)}
\arguments{
  \item{x}{a matrix with HBD probabilities already computes}
  \item{q}{the assumed frequency of the mutation involved in the 
  disease for each individual. Default is 0.00001.}
}
\value{a matrix with FLOD values}
\examples{
require(FEstim)
x <- read.msat.matrix("Example/mapfile.in", "Example/datafile.txt")
x <- festim(x)
x <- HBD.prob(x)
x <- set.FLOD(x)
}
\name{set.HFLOD}
\alias{set.FLOD}
\title{Computes HFLOD probabilities for all the individual's markers}
\description{Returns a matrix with HFLOD computes for each marker of each individual}
\usage{set.HFLOD(x)}
\arguments{
  \item{x}{a matrix with FLOD probabilities already computes (which means it have been passed to festim(), HBD.prob(), set.FLOD() )}
}
\value{a matrix with HFLOD values and alpha max}
\details{
Add the result to the current object
This function uses optimize() to sort the best alpha possible for a given marker
}
\examples{
require(FEstim)
system.file("extdata/map_snp.in", "extdata/data_snp.txt", package="FEstim")
x <- read.msat.matrix("map_snp.in", "data_snp.txt")
x <- festim(x, probs = FALSE)
x <- HBD.prob(x)
x <- set.FLOD(x)
x <- set.HFLOD(x)
}
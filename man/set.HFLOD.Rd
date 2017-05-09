\name{set.HFLOD}
\alias{set.FLOD}
\title{Computes HFLOD probabilities for all the individual's markers}
\description{Returns a matrix with HFLOD computes for each marker of each individual}
\usage{HFLOD.prob(x)}
\arguments{
  \item{x}{a matrix with FLOD probabilities already computes (which means it have been passed to festim(), HBD.prob(), set.FLOD() )}
}
\value{a matrix with HFLOD values and alpha max}
\details{
this function uses optimize() to sort the best alpha possible for a given marker
}
\seealso{ \code{\link{HBD.prob}}, \code{\link{set.FLOD}}, \code{\link{festim}} }
\examples{
require(FEstim)
x <- read.msat.matrix("Example/map_snp.in", "Example/data_snp.txt")
x <- festim(x)
x <- HBD.prob(x)
x <- set.FLOD(x)
x <- set.HFLOD(x)
}
\name{HBD.prob}
\alias{HBD.prob}
\title{Computes HBD probabilities}
\description{Returns a matrix with all the HBD computes for each marker for each individual }
\usage{HBD.prob(x)}
\arguments{
\item{x}{a msat.matrix with a and f computed (which means it have been passed to festim() )}
}
\seealso{\code{\link{festim}}, \code{\link{set.FLOD}}}
\value{A matrix with HBD computes}
\examples{
require(FEstim)
system.file("extdata/in_markall.txt", "extdata/in_datall.txt", package="FEstim")
x <- read.msat.matrix("in_markall.txt", "in_datall.txt")
x <- festim(x)
x <- HBD.prob(x)
}

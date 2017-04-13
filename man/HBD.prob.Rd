\name{HBD.prob}
\title{Computes HBD probabilities for all the individual in the sample}
\description{Returns a matrix with all the HBD computes for each marker for each individual}
\usage{HBD.prob(x)}
\arguments{
\item{x}{a msat.matrix with a and f computed}
}
\value{A matrix with HBD computes}
\examples{
require(FEstim)
x <- read.msat.matrix("Example/in_markallBB.txt", "Example/in_datallBB.txt")
x <- festim(x)
x <- HBD.prob(x)
}

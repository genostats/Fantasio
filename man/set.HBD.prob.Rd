\name{set.HBD.prob}
\alias{set.HBD.prob}
\title{Compute HBD Probabilities}
\description{Add HBD Prob for all the individual in the sample into the object}
\usage{
set.HBD.prob(x)
}
\arguments{ 
	\item{x}{is an object with the coefficient a and f calculated (which means it have been passed to festim() )}
}
\details{before using set.HDB.prob with an object x make sure it have his coefficient a and f estimated with str(nameOftheObject)}
\value{return an object x with a and f estimated}
\seealso{ \code{\link{festim}},  \code{\link{read.msat.matrix}}}
\examples{
require(FEstim)
system.file("extdata/in_markall.txt", "extdata/in_datall.txt", package="FEstim")
x <- read.msat.matrix("Example/in_markall.txt", "Example/in_datall.txt")
x <- festim(x)
x <- HBD.prob(x)
str(x)
}

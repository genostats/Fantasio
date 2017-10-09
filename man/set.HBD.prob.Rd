\name{set.HBD.prob}
\alias{set.HBD.prob}
\title{Compute HBD Probabilities}
\description{Add HBD Prob for all the individual in the sample into the object}
\usage{
set.HBD.prob(x, condition = which(x@p.lrt < 0.05 & x@a <= 1))
}
\arguments{ 
	\item{x}{is an object with the coefficient a and f calculated (which means it have been passed to festim() )}
	\item{condition}{if you want apply a special condition of the data, default is : which(x@p.lrt < 0.05 & x@a <= 1) }
}
\details{before using set.HDB.prob with an object x make sure it have his coefficient a and f estimated with str(nameOftheObject)}
\value{return an object with HBD prob computed}
\examples{
require(FEstim)
system.file("extdata/in_markall.txt", "extdata/in_datall.txt", package="FEstim")
x <- read.msat.matrix("Example/in_markall.txt", "Example/in_datall.txt")
x <- festim(x, probs = FALSE)
x <- set.HBD.prob(x)
str(x)
}

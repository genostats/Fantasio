\name{set.HBD.prob}
\alias{set.HBD.prob}
\title{Add HBD Probabilities To The msat_matrix Object }
\description{Add HBD Prob for all the individual in the sample}
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
x <- set.HBD.prob(x)
str(x)
}

\name{submaps}
\alias{submaps}
\title{Creating several random submaps}
\description{
Return an submap.matrix object that contain multiple submaps and summaries}
\usage{x, n = 100, segment = segments(x), n.cores = 1, epsilon = 1e-3, run.festim = TRUE, proba = TRUE,  verbose = TRUE}
\arguments{
	\item{x}{a bedmatrix}
	\item{n}{the number of submap to create, default is 100}
	\item{segment}{the list of segments outputed by the segments() function}
	\item{n.cores}{the number of core use to create submaps, default is 1}
	\item{run.festim}{if you want to compute a, f, likelihood0, likelihood1, p.lrt, default is TRUE}
	\item{proba}{if you want to compute HBD, FLOD, HFLOD probabilities, default is TRUE}
	\item{verbose}{if you want info from festim() function, default is TRUE }
	}
\details{
If you want to compute HBD, FLOD, HFLOD probabilities make sure that run.festim is TRUE.
}
\value{Return an submap.matrix object with all the submaps created and some summaries.}
\examples{
require(FEstim)
system.file("extdata/geb.bim", "extdata/geb.fam", "extdata/geb.bed", package="FEstim")
x <- read.bed.matrix("geb")
s <- segments(x)
h <- submaps(x, 5, s)  
h
str(h) #to see all the submap's infos
str(h[n]) #with n a particular submap
}

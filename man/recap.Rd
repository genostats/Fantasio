% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recap.r
\name{recap}
\alias{recap}
\title{HBD/FLOD recap}
\usage{
recap(submaps, recap.by.segments = FALSE, list.id)
}
\arguments{
\item{submaps}{a list.submaps object}

\item{recap.by.segments}{: whether you want the recap by segments or snps}

\item{list.id}{: a list individual}
}
\value{
This function returns a list of dataframes.
}
\description{
This function creates HBD and FLOD recap dataframe.
}
\details{
For each individual and each marker the function computes

the mean value of every HBD probabilities computed through the submaps.

The same is done with the FLOD scores.

This function returns a list of two dataframes with HBD and FLOD in it.
}
\examples{
 
bedMatrix <- read.bed.matrix("yourFile")
segmentList <- createSegmentsListByHotspots(bedMatrix)
submaps <- makeSubmapsByHotspots(bedMatrix, 10, segmentList)  
ROH.plot.id(yourROHfile, submaps, FID, IID)
}
\seealso{
Fantasio

makeSubmapsBySnps

createSegmentsListByHotspots

festim

set.HBD.prob

set.FLOD

set.HFLOD

HBD.segments
}
\name{read matrix}
\alias{read.msat.matrix}
\title{create an msat.matrix from mapfile and datafile}
\description{
This function takes in arguments two files : a mapfile and a datafile, 
thus it creates all the arguments necessary for the function msat_matrix()}
\usage{ read.msat.matrix("mapfile.in", "datafile.in")}
\arguments{
	\item{mapfile.in}{a mapfile  sample}
	\item{datafile.in}{a datafile of the sample}
	}
\details{
The function use informations from both files to create an msat.matrix objet which resume all information from previous files.
}
\value{return an msat_matrix object}
\examples{
system.file("extdata/map_msat", "extdata/data_msat.in", package="FEstim")
x <- read.msat.matrix("map_msat.in", "data_msat.in")
x
str(x) #information of your msat.matrix
}

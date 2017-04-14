\name{read matrix}
\alias{read.msat.matrix}
\title{Create All Arguments For the function msat_matrix()}
\description{
This function takes in arguments two files : a mapfile and a datafile, 
thus it creates all the arguments necessary for the function msat_matrix()}
\usage{ read.msat.matrix("mapfile.in", "datafile.in")}
\arguments{
	\item{mapfile.in}{a mapfile  sample}
	\item{datafile.in}{a datafile of the sample}
	}
\details{
the function read.msat.matrix() creates the following arguments:
\itemize{
  \item ncol  = number of genetics marker
  \item nrow  = number of individuals
  \item ped   = a data.frame with the following columns : famId, id, fatherId, motherId, sex,phenotype
  \item msat  = a matrix of genotypes
  \item map   = a data.frame with the following columns : indivId, chromosomes, distance
  \item freq  = a matrix of alleles frequencies
}


After creating them the function read.msat.matrix passes all the variables in arguments for the 
function msat_matrix which creates an object with the information of the sample (use str(nameOftheObject) to watch all the information of your object)
}
\value{return an msat_matrix object}
\seealso{ \code{\link{msat.matrix}}}
\examples{
  x <- read.msat.matrix("map_msat.in", "data_msat.in")
  x
  str(x) #information of your msat.matrix
}

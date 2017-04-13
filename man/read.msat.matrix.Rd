\name{read.msat.matrix}
\title{Create All Arguments For Function msat_matrix()}
\description{
This function takes in arguments two files : a mapfile and a datafile, 
thus it creates all the arguments necessary for the function msat_matrix()}
\usage{ read.msat.matrix("mapfile.in", "datafile.in")}
\arguments{
	\item{mapfile.in}{a mapfile of the sample}
	\item{datafile.in}{a datafile of the sample}
\details{
the function read.msat.matrix() creates the following arguments:
ncol = number of genetics marker
nrow = number of individuals
ped  = a data.frame with the following columns : famId, id, fatherId, motherId, sex, phenotypes
msat = a matrix of genotypes
map  = a data.frame with the following columns : indivId, chromosomes, distance
freq = a matrix of alleles frequencies

After creating them the function read.msat.matrix passes all the variables in arguments for the function msat_matrix
which creates an object with the information of the sample (number of individual and number of markers)
}
\value{return an msat_matrix object}
\seealso{msat_matrix}
\exemples{
x <- read.msat.matrix("map_msat.in", "data_msat.in")
x
str(x)
}

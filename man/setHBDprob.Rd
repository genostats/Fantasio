% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setHBDprob.r
\name{setHBDprob}
\alias{setHBDprob}
\title{Computation of HBD probabilities}
\usage{
setHBDprob(submaps, list.id, quality = 95)
}
\arguments{
\item{submaps}{A atlas object}

\item{list.id}{you can either :
- ignore this parameter if you want to compute HBD, FLOD and HFLOD 
  for individuals who are considerated INBRED and with a QUALITY
  greater or equal to 95%}
- enter a list of individual for a computation of HBD, FLOD score HFLOD score for them
- use "all" for a computation of HBD, FLOD score and HFLOD score for every individual}

\item{quality}{The minimum percentage use to assume if a submap is valid (default is 95)}
}
\value{
the atlas object with each HBD.prob slot of each submaps in the slot submaps_list computed
}
\description{
This function is used to compute HBD probabilities on individuals in a sample
}
\details{
This function iterates over the slots submaps_list of the atlas object.

For each submaps in the slots submaps_list of the object, the slot HBD.prob will be filled with a matrix of dimension : number_inidividual x number_of_markers

By default the function only computes HBD probabilities for INBRED individuals and with a quality equal or greater than 95%. However if you pass the keyword "all" to 
the list.id argument, this function will then computes HBD probabilities for all the individuals in your data. If you want a specific individual, then give a vector containing the family id 
and the individual id separated by an underscore to the list.id argument.
}
\seealso{
setFLOD

setHFLOD
}

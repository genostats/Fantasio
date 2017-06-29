\name{submap}
\alias{submap}
\title{Creating several random submap}
\description{
Return an arraylist that contains multiple submap created thanks to an bedmatrix.}
\usage{submap(x = bedmatrix, n = 100, intensity = 10, hotspot_version = "hg17")}
\arguments{
	\item{x}{a bedmatrix}
	\item{n}{the number of submap to create, default is 100}
	\item{intensity}{the minimal intensity with which the hotspot has to be chosen, default is 10}
	\item{hotspot_version}{the file in which hotspot has to be chosen, default is hg17}
	}
\details{
This function use a file, "hotspot_hg", given by the package that gives hotspots of recombination for the genome, with the "intensity"" for condition the function finds all the hotspots which have an intensity equal or higher to "intensity", with all the hotspots selected the function creates a list of segment.

After that the function finds all the markers between two hotspots (that is to say a segment) and for all the genome.

To put it in a nutshell the function select a random marker in all the segments (delimitate by two hotspots) and creates a submap with all the random markers choosen.
}
\value{Return an arraylist with all the submap created.}
\examples{

system.file("extdata/geb.bim", "extdata/geb.fam", "extdata/geb.bed", package="FEstim")

  x <- read.bed.matrix("geb")
  
  h <- submap(x, 5)
  
  str(h) #to see all the submap's infos
  str(h[n]) #to see the infos of a particular submap 
}

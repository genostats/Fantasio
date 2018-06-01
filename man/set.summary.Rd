\name{setSummary}
\alias{setSummary}
\title{Create all the summary for the submaps.matrix object}
\description{
 Use this function to create all the summary about the submap that you had created.
 }
\usage{ h <- setSummary(h, x)}
\arguments{
  \item{h}{a submap.matrix object}
  \item{x}{a bedmatrix}
}
\details{
  This function allows you to create with one line all the summary for the submap.matrix object
}
\value{return the submap.matrix objet with summary slots completed}
\examples{
  require(FEstim)
  system.file("geb")    #a changer
  x <- read.bed.matrix(geb)
  s <- segments(x)
  h <- submaps(x, 5, s)
  h <- setSummary(h, x)
  str(h)
}

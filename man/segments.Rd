\name{segments}
\alias{segments}
\title{ Create a list of segments }
\description{This function creates a list of segments witch are delimitation of the genome, each segment contains markers. }
\usage{ 
  s <- segments(x)
}

\arguments{
  \item{x}{a bed.matrix}
}

\value{
  An hotspot.segments object.
}

\examples{
  require(FEstim)
  data("extdata/bed/bim/fam.data")
  x <- read.bed.matrix("Exemple/geb")
  s <- segments(x)
  s
}

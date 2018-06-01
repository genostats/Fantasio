\name{segment.summary}
\alias{segment.summary}
\title{ Information from the results of segments() function. }
\description{ create a dataframe with information about the segment that has been created by the segments() function. }
\usage{ 
 df <- segment.summary(s) 
}

\arguments{
  \item{s}{an segment.hotspot object created by the segments() function}
}

\value{
  A dataframe with three columns : 
  \itemize{
    \item{chromosom_number}
    \item{number_of_segments}
    \item{number_of_markers}
  }
}

\examples{
  require(FEstim)
  data("extdata/bed/bim/fam.data")
  x <- read.bed.matrix("Exemple/geb")
  s <- segments(x)
  df <- segment.summary(s)
  df
}

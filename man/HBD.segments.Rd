\name{HBD.segments }
\alias{HBD.segments }
\title{Get HBD segments for the data}
\description{
  Use this function if you want to find HBD segments on your data.
}
\usage{HBD_segments <- HBD.segments(submaps, HBD.recap, n.consecutive.marker = 5, threshold = 0.5, by_segments=F)}
\arguments{
  \item{submaps}{a submap.matrix}
  \item{HBD.recap}{a matrix created by HBD.recap function}
  \item{n.consecutive.marker}{the number of consecutive marker to take and form a segment, default is 5}
  \item{threshold}{the threshold to use when picking a marker, default is 0.5}
  \item{by_segments}{whether to considerate a marker like a whole segment or not, if you used HBD_recap by_segment put by_segment to true, default is false}
}

\value{Return a list of dataframe all the individual with their segments with the following columns :
  \itemize{
    \item{individual}{the individual ID}
    \item{status}{the individual status}
    \item{start}{the beginning of the segment}
    \item{end}{the end of the segment}
    \item{size}{the number of marker in the segment}
    \item{chromosome}{the chromosome which belong the segment}
    \item{start_pos}{the beginning of the segment in cM}
    \item{end_pos}{the end of the segment in cM}
    \item{start_dist}{the beginning of the segment in bp}
    \item{end_dist}{the end of the segment in bp}
  }
}
\examples{
  require(FEstim)
  system.file("geb")    #a changer
  x <- read.bed.matrix(geb)
  s <- segments(x)
  h <- submaps(x, 5, s)
  HBD_recap <- HBD.recap(h)
  HBD_segments <- HBD.segments(h, HBD_recap, threshold = 0.01)
  HBD_segments[1:10, 1:10]
  
  ###if by segment
  
  HBD_recap <- HBD.recap(h, by_segments=T) 
  HBD_segments <- HBD.segments(h, HBD_recap, threshold = 0.01, by_segment=T)
  HBD_segments[1:10, 1:10]
}

\name{HBD.plot.chr}
\alias{HBD.plot.chr}
\title{plot HBD segment for a given chromosome}
\description{
  Use this function if you want to plot HBD segments for a given chromosome.
}
\usage{HBD.plot.id(HBD.recap, HBD.segments, distance="cM", chr, listid="empty", regions="empty", outfile="empty", save_img=F)
\arguments{
  \item{HBD.recap}{a matrix created by HBD.recap function}
  \item{HBD.segments}{a matrix created by HBD.segments function}
  \item{distance}{the unit of the distance between markers, default is cM}
  \item{chr}{the chromosome you want to plot}
  \item{listid}{Allows the user to specify the individuals to select, The file mylist.list should
contain one line per individual with family and individual identifiers Default is empty.}
  \item{regions}{Allows the user to specify regions of interest to be
highlighted in green on the plot. The file myregions.txt should contain one line per region, with its
chromosome and its start and end positions (in cM (default) or bp if the distance is in bp).}
 \item{outfile}{the name of the outfile, default is empty}
 \item{save_img}{whether the plot should be save or not, defautl is false}
}

\value{Return a plot of segment for the given chromosome}
}
\examples{
  require(FEstim)
  system.file("geb")    #a changer
  x <- read.bed.matrix(geb)
  s <- segments(x)
  h <- submaps(x, 5, s)
  HBD_recap <- HBD.recap(h)
  HBD_segments <- HBD.segments(h, HBD_recap, threshold = 0.01)
  HBD.plot.chr(HBD_recap, HBD_segments, chr=1)
}

\name{ROH.plot.id}
\alias{ROH.plot.id}
\title{plot ROH segment for a given individual}
\description{
  Use this function if you want to plot HBD segments for a given individual.
}
\usage{HBD.plot.id(ROHfile, map, distance="cM", regions="empty", outfile="empty", fid, iid, save_file=F  )}
  \arguments{
    \item{ROHfile}{an ROH file}
    \item{map}{a mapfile}
    \item{distance}{the unit of the distance between markers, default is cM}
    \item{chr}{the chromosome you want to plot}
    \item{listid}{Allows the user to specify the individuals to select, The file mylist.list should
      contain one line per individual with family and individual identifiers Default is empty.}
    \item{regions}{Allows the user to specify regions of interest to be
      highlighted in green on the plot. The file myregions.txt should contain one line per region, with its
      chromosome and its start and end positions (in cM (default) or bp if the distance is in bp).}
    \item{outfile}{the name of the outfile, default is empty}
    \item{save_img}{whether the plot should be save or not, defautl is false}
    \item{fid}{family identifier}
    \item{iid}{individual indentifier}
  }
  
  \value{Return a plot of segment for the given individual}

\examples{
}

\name{ROH.plot.chr}
\alias{ROH.plot.chr}
\title{plot ROh segment for a given chromosome}
\description{
  Use this function if you want to plot ROH segments for a given chromosome.
}
\usage{ROH.plot.chr(ROHfile, chr, outfile, listid="empty", regions="empty", distance="cM", map, save_file=F)}
  \arguments{
    \item{ROHfile}{a ROHfile}
    \item{chr}{the chromosome you want the plot from}
    \item{outfile}{the name of the outfile if save_file=T}

    
    \item{listid}{Allows the user to specify the individuals to select. The file mylist.list should
      contain one line per individual with family and individual identifiers id. Default is empty.}
      
    \item{regions}{Allows the user to specify regions of interest to be
      highlighted in green on the plot. The file myregions.txt should contain one line per region, with its
      chromosome and its start and end positions (in cM (default) or bp if the distance is in bp).}
    \item{distance}{the unit distance between marker, default is cM}
    \item{map}{a map file}
    \item{save_img}{whether the plot should be save or not, defautl is false}
  }
  
  \value{Return a plot of segment for the given chromosome}
}
\examples{

}

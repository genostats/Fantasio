##################################################################################
#This function creates HBD segments used for creating plots                      #
#                                                                                #
#!!! submaps : a list object                                                     #                                       
#!!! n.consecutive.marker : the number of consecutive segments/markers with an   #
#    HBD probabiliy >= threshold                                                 #
#!!! threshold : a threshold for n.consecutive.marker                            #
#!!! by_segments : whether we want the segments by snps or segments              #
#                                                                                #
#*** return a list of dataframe with HBD.segment for each individuals            #
##################################################################################

HBD.segments <- function(submaps, n.consecutive.marker = 5, threshold = 0.5, by_segments=FALSE)
{

  if(by_segments)
  {
    HBD.segments.by.segments(submaps, submaps@HBD_recap, n.consecutive.marker, threshold)
  }
  else
  {
    HBD.segments.by.snps(submaps, submaps@HBD_recap, n.consecutive.marker, threshold)
  }
}


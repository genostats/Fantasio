setSummary <- function(h)
{
  h@likelihood_summary <- submap.likelihood(h@atlas)
  h@estimation_summary <- submap.estim(h@atlas)
  h@marker_summary <- summary.marker(h@atlas, h@bedmatrix)
  h@submap_summary <- submap.summary(h@atlas)
  h@HBD_recap <- HBD.recap(h)
  h@HBD_segments <- HBD.segments(h,h@HBD_recap, threshold=0.001) 
  h
}
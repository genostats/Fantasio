setSummary <- function(h, bedmatrix)
{
  h@likelihood_summary <- submap.likelihood(h)
  h@estimation_summary <- submap.estim(h)
  h@marker_summary <- summary.marker(h, bedmatrix)
  h@submap_summary <- submap.summary(h)
  h
}
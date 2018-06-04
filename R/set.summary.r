setSummary <- function(submaps, list.id, run_a_f = TRUE, probs = TRUE, by_segments=FALSE, q=1e-4, threshold=0.5, quality=95, n.consecutive.marker=5)
{
  if(run_a_f)
  {
    submaps@likelihood_summary <- submapLikelihood(submaps@atlas)
  	submaps@estimation_summary <- submapEstim(submaps@atlas)
  	submaps@marker_summary <- summaryMarker(submaps@atlas, submaps@bedmatrix)
  	submaps@submap_summary <- submapSummary(submaps@atlas)
  }
  
  if(run_a_f && probs)
  {
    l1 <- set.HBD.prob(submaps, list.id=list.id, quality=quality)
    submaps <- l1[[1]]
    submaps <- set.FLOD(submaps=submaps, condition=l1[[2]], q=q)
    
    if(!(class(submaps@atlas[[1]])[1] == "snps.matrix" & by_segments))
    {
      l2 <- recap(submaps, by_segments=by_segments, list.id=list.id)
  	  submaps@HBD_recap <- l2[[1]]
  	  submaps@FLOD_recap <- l2[[2]]
  	  submaps@HBD_segments <- HBD.segments(submaps, threshold=threshold, n.consecutive.marker=n.consecutive.marker) 
  	  submaps@HFLOD <- set.HFLOD(submaps)
    }
  	
  }  
  submaps
}
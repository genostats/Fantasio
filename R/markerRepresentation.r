#' Markers selected through the submaps
#' 
#' This function gives for each marker 
#' the number of times it has been selected in a submap.
#' 
#' @param atlas an atlas
#' 
#' @return A data frame with columns: 
#' \describe{
#'  \item{id}{id of the snps}
#'  \item{chr}{chromosome}
#'  \item{pos}{position in bp}
#'  \item{dist}{position in cM}
#'  \item{Freq}{number of times the marker has been picked (0 for markers never selected in any submap)}
#' }
#' 
#' @seealso setSummary
#' 
#' @export
markerRepresentation <- function(atlas) {
    if (class(atlas)[1] != "atlas")
        stop("need an atlas")
    res <- data.frame(id = atlas@bedmatrix@snps$id, chr = atlas@bedmatrix@snps$chr, pos = atlas@bedmatrix@snps$pos, 
                 dist = atlas@bedmatrix@snps$dist, Freq = 0, stringsAsFactors = FALSE)

    submaps <- atlas@submaps_list
    id <- unlist(sapply(submaps, function(x) x@map$id, simplify = FALSE))
    b <- as.data.frame(table(id), stringsAsFactors = FALSE)
    m <- match( b$id, res$id )
    res$Freq[m] <- b$Freq
    res
}


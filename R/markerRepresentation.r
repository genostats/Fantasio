#' Markers selected through the submaps
#' 
#' This function gives for each marker (picked at least once) 
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
#'  \item{Freq}{number of times the marker has been picked}
#' }
#' 
#' @seealso setSummary
#' 
#' @export
markerRepresentation <- function(atlas) {
    submaps <- atlas@submaps_list
    if (class(atlas)[1] != "atlas")
        stop("need an atlas")
    id <- unlist(sapply(submaps, function(x) x@map$id, simplify = FALSE))
    b <- as.data.frame(table(id), stringsAsFactors = FALSE)
    m <- match( b$id, atlas@bedmatrix@snps$id )
    b$chr <- atlas@bedmatrix@snps$chr[m]
    b$pos <- atlas@bedmatrix@snps$pos[m]
    b$dist <- atlas@bedmatrix@snps$dist[m]
    b <- b[ order(b$chr, b$pos), ]
    b <- b[, c(1,3,4,5,2)]
    rownames(b) <- NULL
    b
}


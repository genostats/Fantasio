setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClass("msat.matrix", representation(
           # input slots
           ncol = 'numeric',            # nb msats
           nrow = 'numeric',            # nb inds
            ped = 'data.frame',         # le contenu typique d'un .fam / 6 premiC(res cols d'un .ped
           msat = 'matrix',             # matrice de gC)notypes (nb msats x 2 nb inds)
            map = 'data.frame',         # id, chr, distance
           freq = 'matrix',             # matrice de frC)quences allC)liques (nb msat x max allele)
        epsilon = 'numericOrNULL',      # la valeur d'epsilon utilisC)e pour calculer la prC)cC)dente
           # output slots
       log.freq = 'matrix',             # log de la prC)cC)dente (il est utile qu'elle soit prC)calculC)e)
     delta.dist = 'numericOrNULL',      # diff(map$distance) + faire attention au chgt de chr
      log.emiss = 'matrixOrNULL',       # matrice de logs proba d'emission si statut = 0 ou 1 (2 nb inds x nb msats) (mC*me commentaire)
              a = 'numeric',            # valeurs de a et f estimC)es par la fonction festim
              f = 'numeric',
    likelihood0 = 'numeric',            # vraisemblance sous H0 (f = 0)
    likelihood1 = 'numeric',            # vraisemblance sous H1, ie aux valeurs calculés 
          p.lrt = 'numeric',            # le likelihood ratio test
       HBD.prob = 'matrix',             # proba d'etre HBD = 1 ; un individu par colonne : dim = (nb inds x nb msats)
           FLOD = 'matrix',             # matrice des FLOD scores  dim = (nb inds x nb msats)
          HFLOD = 'matrix',              # matrice des HFLOD scores dim = (nb msats x 2)
         Submap = 'matrix'
) )

setMethod('initialize', signature='msat.matrix', definition=function(.Object,ncol,nrow,ped,msat,map,freq) {
  #if(nrow(freq) != ncol) stop("freq dimension mismatch")
  #if(nrow(msat) != ncol | ncol(msat) != 2*nrow) stop("msat dimension mismatch")
  if(nrow(ped) != nrow) stop("ped dimension mismatch")
  w <- is.na(msat)
  if(any(w)) {
    warn('NA msats replaced by zeros')
    msat[w] <- 0
  }
  if(max(msat) > ncol(freq)) stop("allele exceeding number of columns of freq")

  ## TODO 
  # vC)rifier que dans map, les marqueurs sont bien ordonnC)s sur chaque chromosome...
  # et que les chromosomes ne sont pas en vrac

  # mC*me code dans random.msat
  delta.dist <- diff(map$distance)
  # traiter proprement les changements de chromosomes
  I <- cumsum(rle(map$chr)$length)
  I <- I[-length(I)]
  delta.dist[I] <- -1 

  .Object@ncol     <- ncol
  .Object@nrow     <- nrow
  .Object@ped      <- ped
  .Object@msat     <- msat
  .Object@map      <- map 
  .Object@freq     <- freq
  .Object@log.freq <- log(freq)
  .Object@delta.dist <- delta.dist
  .Object
})

# x = une msat.matrix
set.log.emiss <- function(x, epsilon = 1e-3) {
  if(is.null(x@log.freq)) x@log.freq <- log(x@freq)
  logEmiss <- matrix(ncol = x@ncol, nrow = 2*x@nrow)
  for(i in 1:x@nrow) {
    Y1 <- x@msat[,2*i-1]
    Y2 <- x@msat[,2*i]
    logEmiss[ c(2*i-1,2*i) ,] <- .Call('festim_logEmiss', PACKAGE = "FEstim", Y1, Y2, x@log.freq, epsilon)
  }
  x@log.emiss <- logEmiss
  x@epsilon <- epsilon
  x
}


setGeneric('dim')
setMethod("dim", signature = "msat.matrix", function(x) c(x@nrow, x@ncol) )

setMethod('show', signature("msat.matrix"),
  function(object) {
    cat('A msat.matrix with ', nrow(object), ' individual(s) and ', ncol(object), " markers\n")
  }
)


setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClass("f.matrix", representation(
          # input slots
          ncol = 'numeric',              # nb loci
          nrow = 'numeric',              # nb inds
          ped = 'data.frame',            # le contenu typique d'un .fam / 6 premieres cols d'un .ped
          map = 'data.frame',            # id, chr, distance
          epsilon = 'numeric',     # la valeur d'epsilon utilisee pour calculer les log emission
          # dérive de map, précalculé à l'initialisation
          delta.dist = 'numericOrNULL',  # diff(map$distance) + faire attention au chgt de chr
          # créés à partir d'une msat matrix / d'une bed matrix...
          log.emiss = 'matrixOrNULL',    # matrice de logs proba d'emission si statut = 0 ou 1 (2 nb inds x nb msats) (meme commentaire)
          # output slots
          a = 'numeric',                 # valeurs de a et f estimees par la fonction festim
          f = 'numeric',
          likelihood0 = 'numeric',       # vraisemblance sous H0 (f = 0)
          likelihood1 = 'numeric',       # vraisemblance sous H1, ie aux valeurs calculés 
          p.lrt = 'numeric',             # le likelihood ratio test
          HBD.prob = 'matrix',           # proba d'etre HBD = 1 ; un individu par colonne : dim = (nb inds x nb msats)
          FLOD = 'matrix',               # matrice des FLOD scores  dim = (nb inds x nb msats)
          HFLOD = 'matrix'               # matrice des HFLOD scores dim = (nb msats x 2)
))

setMethod('initialize', signature='f.matrix', definition=function(.Object, ncol, nrow, ped, map, log.emiss, epsilon=1e-3) {
  if(nrow(ped) != nrow) stop("ped dimension mismatch")
  if(nrow(map) != ncol) stop("map dimension mismatch")
 
  ## TODO 
  # verifier que dans map, les marqueurs sont bien ordonnes sur chaque chromosome...
  # et que les chromosomes ne sont pas en vrac
 
  # meme code dans random.msat
  delta.dist <- diff(map$distance)
  # traiter proprement les changements de chromosomes
  I <- cumsum(rle(map$chr)$length)
  I <- I[-length(I)]
  delta.dist[I] <- -1 
  
  .Object@ncol     <- ncol
  .Object@nrow     <- nrow
  .Object@ped      <- ped
  .Object@map      <- map 
  .Object@delta.dist <- delta.dist
  .Object@log.emiss <- log.emiss
  .Object@epsilon   <- epsilon
  .Object
})

setGeneric('dim')
setMethod("dim", signature = "f.matrix", function(x) c(x@nrow, x@ncol) )

setMethod('show', signature("f.matrix"),
          function(object) {
            cat('A f.matrix with ', nrow(object), ' individual(s) and ', ncol(object), " markers\n")
          }
)

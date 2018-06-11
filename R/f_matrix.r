##################################################################################
#This is the mother class for any object of this package                         #
#                                                                                #
##################################################################################


setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClass("f.matrix", representation(
          # input slots
          ncol = 'numeric',              # number of loci
          nrow = 'numeric',              # number of individual
          ped = 'data.frame',            # the first 6 columns of a .ped
          map = 'data.frame',            # id, chr, distance
          epsilon = 'numeric',           # value of epsilon = genotyping error rare, to use for computing log emission precalculated at initialisation
          delta.dist = 'numericOrNULL',  # diff(map$distance) + faire attention au chgt de chr
          
          # slot created from an msatmatrix/bed.matrix
          
          log.emiss = 'matrixOrNULL',    # matrix of log proba emission if statut = 0 or 1(2 nb inds x nb msats) 
          
          # output slots
          
          a = 'numeric',                 # value of a and f estimated by festim
          f = 'numeric',
          likelihood0 = 'numeric',       # likelihood under H0 (f = 0)
          likelihood1 = 'numeric',       # likelihood under H1 
          p.lrt = 'numeric',             # likelihood ratio test
          HBD.prob = 'matrix',           # proba HBD = 1 ; one individual per column : dim = (nb inds x nb msats)
          FLOD = 'matrix'               # matrix of FLOD scores dim = (nb inds x nb msats)
          #HFLOD = 'matrix'               # matrix of HFLOD scores dim = (nb msats x 2)
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

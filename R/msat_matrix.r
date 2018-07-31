##################################################################################
#This function  is uses to create a msat.matrix object similar to a submap.matrix#
#object                                                                          #
#                                                                                #
##################################################################################


setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClass("msat.matrix", representation(
           # input slots
           msat = 'matrix',             # genotype matrix (nb msats x 2 nb inds)
           freq = 'matrix',             # allelic frequencies matrix (nb msat x max allele)
           # output slots
       log.freq = 'matrix'              # previous log 
), contains = "f.matrix" )

setMethod('initialize', signature='msat.matrix', definition=function(.Object, ncol, nrow, ped, msat, map, freq, epsilon = 1e-3) {
  if(nrow(freq) != ncol) stop("freq dimension mismatch")
  if(nrow(msat) != ncol | ncol(msat) != 2*nrow) stop("msat dimension mismatch")
  w <- is.na(msat)
  if(any(w)) {
    warning('NA msats replaced by zeros')
    msat[w] <- 0
  }
  if(max(msat) > ncol(freq)) stop("allele exceeding number of columns of freq")



  .Object@msat     <- msat
  .Object@freq     <- freq
  .Object@log.freq <- log(freq)
  .Object@ncol     <- ncol  
  .Object@nrow     <- nrow
  .Object@epsilon  <- epsilon
  log.emiss <- msat.log.emiss(.Object)
  callNextMethod(.Object, ncol, nrow, ped, map, log.emiss, epsilon)
})

# x = une msat.matrix
msat.log.emiss <- function(x) {
  if(is.null(x@log.freq)) x@log.freq <- log(x@freq)
  logEmiss <- matrix(ncol = x@ncol, nrow = 2*x@nrow)
  for(i in 1:x@nrow) {
    Y1 <- x@msat[,2*i-1]
    Y2 <- x@msat[,2*i]
    logEmiss[ c(2*i-1,2*i) ,] <- .Call('festim_logEmiss', PACKAGE = "Fantasio", Y1, Y2, x@log.freq, x@epsilon)
  }
  logEmiss
}

setMethod('show', signature("msat.matrix"),
  function(object) {
    cat('A msat.matrix with ', nrow(object), ' individual(s) and ', ncol(object), " markers\n")
  }
)


set.HFLOD <- function(x, step = 0.01)
{
  #initialiser alpha
 alpha <- matrix(0.0)
 nb <- 0
 i <- 0
 while(nb * step < 1)
 {
   alpha[i] <- nb * step
   nb <- nb + 1
   i <- i + 1
 }
 alpha[i] <- 1
  

 
  h <- function(alpha)
  {
      sum <- c()
     
       for (i in 1:x@nrow)
       {
          if(x@f[i]==0) next
          v <- sum(log10(alpha*exp(x@FLOD[i,] *log(10))+(1-alpha)))
          #v <- sum(log(alpha * (exp(x@FLOD[i,]*log(10))-1) + 1)/log(10))
          sum <- c(sum,v)
       }
    return(sum)
  }
}
# Error in optim(par = c(0, 1), fn = h) : 
#objective function in optim evaluates to length 8 not 1
#besoin de retourner un vecteur de longueur 1

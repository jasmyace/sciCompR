eye <- function(n){
  M <- matrix(rep(0,n*n),nrow=n,ncol=n)
  for(i in 1:n){
    M[i,i] <- 1 
  }
  return(M)
}
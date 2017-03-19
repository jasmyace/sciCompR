diagM <- function(A){
  D <- matrix(rep(0,nrow(A)*nrow(A)),nrow=nrow(A),ncol=ncol(A))
  for(i in 1:nrow(A)){
    D[i,i] <- A[i,i]
  }
  return(D)
}

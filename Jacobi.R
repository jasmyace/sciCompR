
#   ---- Jacobi.

x0 <- matrix(c(0,0,0),nrow=3,ncol=1)
A <- matrix(c(10,1,-1,1,15,1,-1,1,20),nrow=3,ncol=3,byrow=TRUE)
b <- c(18,-12,17)

eye <- function(n){
  M <- matrix(rep(0,n*n),nrow=n,ncol=n)
  for(i in 1:n){
    M[i,i] <- 1 
  }
  return(M)
}

diagM <- function(A){
  D <- matrix(rep(0,nrow(A)*nrow(A)),nrow=nrow(A),ncol=ncol(A))
  for(i in 1:nrow(A)){
    D[i,i] <- A[i,i]
  }
  return(D)
}
diagM(A)



jacobi <- function(A,b,method="m",epsilon=1e-3,maxIters=1000){
  n <- nrow(A)
  d <- diag(A)
  convergence <- FALSE
  
  for(j in 1:n){
    x[j] <- 0
  }
  
  i <- 1
  
  repeat{
    
    if( (norm < epsilon & i > 2) | i == maxIters ){
      convergence <- TRUE
      return(list(x=x,iters=i - 1,norm=norm))
      break
    } else {

      if(method == "v"){

        xh <- x
        for(j in 1:n){
          s <- 0
          for(k in 1:n){
            if(k != j){
              s <- s - A[j,k]*xh[k] 
            }
          }
          x[j] <- (s + b[j]) /  d[j]
        }
      } else if(method == "m"){
        D <- diagM(A)
        if(i == 1){
          for(j in 1:n){
            x[j] <- 0
          }
        } else {
          xh <- x
          x <- solve(D,(D - A) %*% x + b)
        } 
      }
      norm <- norm(x - xh,"f")
      cat(paste0("Iteration ",i,":  Norm = ",round(norm,3),", Epsilon = ",epsilon,", Vector = ",t(x),"\n"))
      i <- i + 1
    }
  }
}





gs <- function(A,b,method="v",epsilon=1e-3,maxIters=1000){
  n <- nrow(A)
  d <- diag(A)
  convergence <- FALSE
  
  for(j in 1:n){
    x[j] <- 0
  }
  
  i <- 1
  
  repeat{
    
    if( (norm < epsilon & i > 2) | i == maxIters ){
      convergence <- TRUE
      return(list(x=x,iters=i - 1,norm=norm))
      break
    } else {
      
      if(method == "v"){
        
        xh <- x
        for(j in 1:n){
          s <- 0
          for(k in 1:n){
            if(k != j){
              s <- s - A[j,k]*x[k] 
            }
          }
          x[j] <- (s + b[j]) /  d[j]
        }
      } else if(method == "m"){
#         D <- diagM(A)
#         if(i == 1){
#           for(j in 1:n){
#             x[j] <- 0
#           }
#         } else {
#           xh <- x
#           x <- solve(D,(D - A) %*% x + b)
#         } 
      } 
      norm <- norm(x - xh,"f")
      cat(paste0("Iteration ",i,":  Norm = ",round(norm,3),", Epsilon = ",epsilon,", Vector = ",t(x),"\n"))
      i <- i + 1
    }
  }
}




A <- matrix(c(3,1,1,1,-5,2,2,1,5),nrow=3,ncol=3,byrow=TRUE)
b <- c(1,-7,10)

jacobi(A,b,method="v",epsilon=10e-6)
gs(A,b,method="v",epsilon=10e-6)



require(matrixcalc)
LU <- lu.decomposition(A)
L <- LU[[1]]
U <- LU[[2]]


A <- matrix(c(3,1,1,1,-5,2,2,1,5),nrow=3,ncol=3,byrow=TRUE)
b <- c(0,5,10)
LUDecomp <- function(A,b,M=20){
  LU <- lu.decomposition(A)
  L <- LU[[1]]
  U <- LU[[2]]
  X1 <- matrix(rep(NA,length(b)*M),nrow=20,ncol=length(b))
  X2 <- matrix(rep(NA,length(b)*M),nrow=20,ncol=length(b))
  for (i in 1:M){
    b[1] <- i
    X1[i,] <- matrix(solve(U,solve(L,b)),ncol=1,nrow=length(b))
    X2[i,] <- matrix(solve(L,solve(U,b)),ncol=1,nrow=length(b))
  }
  return(list(X1=X1,X2=X2))
}
test <- LUDecomp(A,b)




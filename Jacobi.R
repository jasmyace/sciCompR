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
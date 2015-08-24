working.comp <- function(fit){
   
  G <- -(fit$gradient - fit$ps$S.h2)
  H <- fit$hessian - fit$ps$S.h
  
  W.eig <- eigen(H, symmetric=TRUE)
  
  if(min(W.eig$values) < sqrt(.Machine$double.eps)) { pep <- which(W.eig$values < sqrt(.Machine$double.eps)); W.eig$values[pep] <- 0.0000001 }  
  
  srev    <- sqrt(W.eig$val)
  c.W     <- W.eig$vec%*%tcrossprod(diag(srev)  ,W.eig$vec) 
  W.invsr <- W.eig$vec%*%tcrossprod(diag(1/srev),W.eig$vec)
  
  X <- c.W 
  Z <- W.invsr%*%G + X%*%fit$argument 
  
  
  
  
        #op <- options(warn = -1)
        #R <- chol(H, pivot = TRUE)
        #options(op)
        #p <- dim(H)[2]
        #ipiv <- piv <- attr(R, "pivot")
        #ipiv[piv] <- 1:p
        #rank <- attr(R, "rank")
        #ind <- 1:rank
        #if (rank < p) R[(rank + 1):p, ] <- 0
        #R <- X <- R[ipiv, ipiv]
        #f <- X%*%Z
        #Z <- f <- forwardsolve(R, f)
  
  
  
  
          
  rm(G, H, W.eig, srev, c.W, W.invsr)        
          
 list( X = X , Z = Z )

}



























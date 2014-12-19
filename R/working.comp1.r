working.comp1 <- function(fit){
   
  nz <- 0.000000001
  
  G <- -(fit$gradient - fit$ps$S.h2)
  H <- fit$hessian - fit$ps$S.h
  
  W.eig <- eigen(H, symmetric=TRUE)
  
  if(min(W.eig$values) < .Machine$double.eps) { pep <- which(W.eig$values < .Machine$double.eps); W.eig$values[pep] <- nz }  
  
  srev    <- sqrt(W.eig$val)
  c.W     <- W.eig$vec%*%tcrossprod(diag(srev)  ,W.eig$vec) 
  W.invsr <- W.eig$vec%*%tcrossprod(diag(1/srev),W.eig$vec)
  
  X <- c.W 
  Z <- W.invsr%*%G + X%*%fit$argument 
                 
 list( X = X , Z = Z )

}



























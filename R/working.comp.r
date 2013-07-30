working.comp <- function(x,X1=X1,X2=X2,X1.d2=X1.d2,X2.d2=X2.d2,n=n){

  e.par <- x$argument
  e.par <- e.par[-length(e.par)]
  X <- rW.X <- matrix(0,2*n,(X1.d2+X2.d2))
  D <- rW.Z <- matrix(0,2*n,1)
  j <- 1

    for(i in seq(1,(2*n-1),by=2)) {
      X[i,1:X1.d2]                   <- X1[j,]
      X[i+1,(X1.d2+1):(X1.d2+X2.d2)] <- X2[j,]

      D[i,1]   <- x$dl.dbe1[j]
      D[i+1,1] <- x$dl.dbe2[j]

      W <- matrix(c( x$d2l.be1.be1[j],x$d2l.be1.be2[j],     
                     x$d2l.be1.be2[j],x$d2l.be2.be2[j]), 2 , 2 ) 

      W.eig <- eigen(W,symmetric=TRUE)
      c.W   <- W.eig$vec%*%tcrossprod(diag(sqrt(pmax(W.eig$val,.Machine$double.eps))),W.eig$vec) 
      W.inv <- W.eig$vec%*%tcrossprod(diag(1/pmax(W.eig$val,.Machine$double.eps)),    W.eig$vec) 
      rW.X[i:(i+1),]  <- c.W%*%X[i:(i+1),]
      rW.Z[i:(i+1),1] <- c.W%*%( X[i:(i+1),]%*%e.par + W.inv%*%D[i:(i+1),1] )

      j <- j + 1
    }
 list( rW.X=rW.X , rW.Z=rW.Z )
}





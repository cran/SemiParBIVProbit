working.comp <- function(x,X1=X1,X2=X2,X1.d2=X1.d2,X2.d2=X2.d2,n=n){
  #x=fit
  der.1 <- x$dl.dbe1
  der.2 <- x$dl.dbe2
  der.3 <- x$dl.drho
  e.par <- x$argument
  W <- c.W <- i.W <- matrix(0,3,3)
  X <- rW.X <- matrix(0,3*n,(X1.d2+X2.d2+1))
  D <- iW.D <- Z <- rW.Z <- matrix(0,3*n,1)
  j <- 1
  #n.e=3
  #i=1
    for(i in seq(1,(3*n-2),by=3)) {
      X[i,1:X1.d2]                   <- X1[j,]
      X[i+1,(X1.d2+1):(X1.d2+X2.d2)] <- X2[j,]
      X[i+2,(X1.d2+X2.d2+1)]         <- 1
      D[i,1]   <- der.1[j]
      D[i+1,1] <- der.2[j]
      D[i+2,1] <- der.3[j]
      W <- matrix(c( x$d2l.be1.be1[j],x$d2l.be1.be2[j],x$d2l.be1.rho[j],     
                     x$d2l.be1.be2[j],x$d2l.be2.be2[j],x$d2l.be2.rho[j],  
                     x$d2l.be1.rho[j],x$d2l.be2.rho[j],x$d2l.rho.rho[j] ) , 3 , 3 ) 
      #c.W <- chol(W),pivot=TRUE)
      #piv <- order(attr(c.W,"pivot")) 
      #c.W <- c.W[,piv]
      W.eig <- eigen(W)
      c.W   <- W.eig$vec%*%diag(sqrt(pmax(W.eig$val,.Machine$double.eps^0.6)))%*%t(W.eig$vec) 
      W.inv <- W.eig$vec%*%diag(1/pmax(W.eig$val,.Machine$double.eps^0.6))%*%t(W.eig$vec) 
      rW.X[i:(i+2),1:(X1.d2+X2.d2+1)] <- c.W%*%X[i:(i+2),1:(X1.d2+X2.d2+1)]
      rW.Z[i:(i+2),1] <- c.W%*%( X[i:(i+2),]%*%e.par + W.inv%*%D[i:(i+2),1] )
      #solve(W,D[i:(i+2),1]) 
      #print(j)
      j <- j + 1
    }
 list( rW.X=rW.X , rW.Z=rW.Z )
}



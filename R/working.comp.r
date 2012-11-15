working.comp <- function(x,X1=X1,X2=X2,X1.d2=X1.d2,X2.d2=X2.d2,n=n){

  e.par <- x$argument
  X <- rW.X <- matrix(0,3*n,(X1.d2+X2.d2+1))
  D <- rW.Z <- matrix(0,3*n,1)
  j <- 1

    for(i in seq(1,(3*n-2),by=3)) {
      X[i,1:X1.d2]                   <- X1[j,]
      X[i+1,(X1.d2+1):(X1.d2+X2.d2)] <- X2[j,]
      X[i+2,(X1.d2+X2.d2+1)]         <- 1
      D[i,1]   <- x$dl.dbe1[j]
      D[i+1,1] <- x$dl.dbe2[j]
      D[i+2,1] <- x$dl.drho[j]
      W <- matrix(c( x$d2l.be1.be1[j],x$d2l.be1.be2[j],x$d2l.be1.rho[j],     
                     x$d2l.be1.be2[j],x$d2l.be2.be2[j],x$d2l.be2.rho[j],  
                     x$d2l.be1.rho[j],x$d2l.be2.rho[j],x$d2l.rho.rho[j] ) , 3 , 3 ) 

      W.eig <- eigen(W)
      c.W   <- W.eig$vec%*%tcrossprod(diag(sqrt(pmax(W.eig$val,sqrt(.Machine$double.eps)))),W.eig$vec) 
      W.inv <- W.eig$vec%*%tcrossprod(diag(1/pmax(W.eig$val,sqrt(.Machine$double.eps))),W.eig$vec) 
      rW.X[i:(i+2),1:(X1.d2+X2.d2+1)] <- c.W%*%X[i:(i+2),1:(X1.d2+X2.d2+1)]
      rW.Z[i:(i+2),1] <- c.W%*%( X[i:(i+2),]%*%e.par + W.inv%*%D[i:(i+2),1] )

      j <- j + 1
    }
 list( rW.X=rW.X , rW.Z=rW.Z )
}



#working.comp <- function(x,X1=X1,X2=X2,X1.d2=X1.d2,X2.d2=X2.d2,n=n){
#
#  e.par <- x$argument
#  c.W <- W.inv <- list()
#  X <- rW.X <- matrix(0,3*n,(X1.d2+X2.d2+1))
#  D <- iW.D <- Z <- rW.Z <- matrix(0,3*n,1)
#  j <- 1
#
#    for(i in seq(1,(3*n-2),by=3)) {
#      X[i,1:X1.d2]                   <- X1[j,]
#      X[i+1,(X1.d2+1):(X1.d2+X2.d2)] <- X2[j,]
#      X[i+2,(X1.d2+X2.d2+1)]         <- 1
#      D[i,1]   <- x$dl.dbe1[j]
#      D[i+1,1] <- x$dl.dbe2[j]
#      D[i+2,1] <- x$dl.drho[j]
#      W <- matrix(c( x$d2l.be1.be1[j],x$d2l.be1.be2[j],x$d2l.be1.rho[j],     
#                     x$d2l.be1.be2[j],x$d2l.be2.be2[j],x$d2l.be2.rho[j],  
#                     x$d2l.be1.rho[j],x$d2l.be2.rho[j],x$d2l.rho.rho[j] ) , 3 , 3 ) 
#
#      W.eig <- eigen(W)
#      c.W[[j]]   <- W.eig$vec%*%tcrossprod(diag(sqrt(pmax(W.eig$val,sqrt(.Machine$double.eps)))),W.eig$vec) 
#      W.inv[[j]] <- W.eig$vec%*%tcrossprod(diag(1/pmax(W.eig$val,sqrt(.Machine$double.eps))),W.eig$vec) 
#
#      j <- j + 1
#    }
# rW.X <- bdiag(c.W)%*%X
# rW.Z <- bdiag(c.W)%*%( X%*%e.par + bdiag(W.inv)%*%D )
#
# list( rW.X=rW.X , rW.Z=rW.Z )
#}


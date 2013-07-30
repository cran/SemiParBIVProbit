working.compNP <- function(x,X1=X1,X2=X2,X1.d2=X1.d2,X2.d2=X2.d2,n=n,K=NULL){

  X <- rW.X <- matrix(0,2*n,(X1.d2+X2.d2+2*K))
  D <- rW.Z <- matrix(0,2*n,1)
  e.par     <- x$argument
  e.par <- e.par[-length(e.par)]
  j         <- 1

  for(i in seq(1,(2*n-1),by=2)) {  

  	for(u in 1:K){

 	 X[i,1:(X1.d2+K)]                     <- X1(u)[j,]
  	 X[i+1,(X1.d2+1+K):(X1.d2+X2.d2+2*K)] <- X2(u)[j,]
  	 D[i,1]   <- x$dl.dbe1[j,u]
  	 D[i+1,1] <- x$dl.dbe2[j,u]

  	 W <- matrix(c( x$d2l.be1.be1[j,u],x$d2l.be1.be2[j,u],     
             	        x$d2l.be1.be2[j,u],x$d2l.be2.be2[j,u] ) , 2 , 2 ) 

  	 W.eig <- eigen(W,symmetric=TRUE)
  	 c.W   <- W.eig$vec%*%tcrossprod(diag(sqrt(pmax(W.eig$val,sqrt(.Machine$double.eps)))),W.eig$vec) 
   	 W.inv <- W.eig$vec%*%tcrossprod(diag(1/pmax(W.eig$val,sqrt(.Machine$double.eps))),    W.eig$vec) 
      
  	 rW.X[i:(i+1),]  <- rW.X[i:(i+1),]  + x$We[j,u]*c.W%*%X[i:(i+1),]
  	 rW.Z[i:(i+1),1] <- rW.Z[i:(i+1),1] + x$We[j,u]*( c.W%*%(X[i:(i+1),]%*%e.par + W.inv%*%D[i:(i+1),1]) )

  	}

  j <- j + 1

  }

 list( rW.X=rW.X , rW.Z=rW.Z )
}



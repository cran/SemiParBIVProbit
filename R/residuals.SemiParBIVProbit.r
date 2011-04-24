residuals.SemiParBIVProbit <- function(object,...){
  
  der.1 <- object$fit$dl.dbe1
  der.2 <- object$fit$dl.dbe2
  der.3 <- object$fit$dl.drho

  D <- matrix(0,3*object$n,1)
  r.w <- r.p <- matrix(0,object$n,3) 

  j <- 1

    for(i in seq(1,(3*object$n-2),by=3)) {

      D[i,1]   <- der.1[j]
      D[i+1,1] <- der.2[j]
      D[i+2,1] <- der.3[j]

      W <- matrix(c( object$fit$d2l.be1.be1[j],object$fit$d2l.be1.be2[j],object$fit$d2l.be1.rho[j],     
                     object$fit$d2l.be1.be2[j],object$fit$d2l.be2.be2[j],object$fit$d2l.be2.rho[j],  
                     object$fit$d2l.be1.rho[j],object$fit$d2l.be2.rho[j],object$fit$d2l.rho.rho[j] ) , 3 , 3 ) 

      W.eig <- eigen(W)
      W.ins <- W.eig$vec%*%diag(1/sqrt(pmax(W.eig$val,.Machine$double.eps^0.6)))%*%t(W.eig$vec) 
      W.inv <- W.eig$vec%*%diag(1/pmax(W.eig$val,.Machine$double.eps^0.6))%*%t(W.eig$vec) 

      r.p[j,] <- W.ins%*%D[i:(i+2),1]
      r.w[j,] <- W.inv%*%D[i:(i+2),1]
  
      j <- j + 1

    }

 list( r.p=r.p, r.w=r.w )

}
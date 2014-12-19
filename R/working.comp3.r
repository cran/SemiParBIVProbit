working.comp3 <- function(x, VC = VC, myf = myf,  extra.l =  extra.l){

  good  <- x$good
  n     <- sum(as.numeric(good==TRUE))
  sqq   <- seq(1, (3*n-(3-1)), by = 3)
  e.par <- x$argument

  D             <- extra.l$D[good,]
  D[sqq,     1] <- x$dl.dbe1
  D[sqq + 1, 1] <- x$dl.dbe2
  D[sqq + 2, 1] <- x$dl.drho
  
  X                                            <- extra.l$X[good,]
  X[sqq,      1:VC$X1.d2]                      <- VC$X1[good,] 
  X[sqq + 1, (VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)] <- VC$X2[good,]
  X[sqq + 2, (VC$X1.d2 + VC$X2.d2 + 1)]        <- 1
  

  be1be2 <- cbind(x$d2l.be1.be1, x$d2l.be1.be2, x$d2l.be2.be2, x$d2l.be1.rho, x$d2l.be2.rho, x$d2l.rho.rho)

  L.W <- unlist( apply(be1be2, 1, myf ) , recursive = FALSE) 

  W.invsr <- c.W <- list()
  nz <- 0.000000001

    for(i in 1:n){

      s.ch <- sum(as.numeric(as.vector(L.W[[i]])==0))
      
     if(VC$Model=="BSS" && s.ch==8){ 
      
                          if(L.W[[i]][1,1] < .Machine$double.eps) L.W[[i]][1,1] <- nz
                          c.W[[i]] <- W.invsr[[i]] <- matrix(0,3,3)
                          c.W[[i]][1,1] <- sqrt(L.W[[i]][1,1]) 
                          W.invsr[[i]][1,1] <- 1/c.W[[i]][1,1]
                          
                                   }else{
                                   
      W.eig <- eigen(L.W[[i]], symmetric=TRUE)                               
      if(min(W.eig$values) < .Machine$double.eps) { pep <- which(W.eig$values < .Machine$double.eps)
      	                                            W.eig$values[pep] <- nz  
      	                                           }                                   
                                        }
      
      if( !(VC$Model=="BSS" && s.ch==8) ){
      
      srev <- sqrt(W.eig$val)
      c.W[[i]]     <- W.eig$vec%*%tcrossprod(diag(srev)  ,W.eig$vec) 
      W.invsr[[i]] <- W.eig$vec%*%tcrossprod(diag(1/srev),W.eig$vec)
      
      }

                  }

      rW.X <- as.matrix( bdiag(c.W) %*% X ) 
      rW.Z <- as.matrix( rW.X%*%e.par + bdiag(W.invsr) %*% D ) 

rm(e.par, c.W, W.invsr, D, be1be2, L.W )
if(VC$gc.l == TRUE) gc()

                   
 list( rW.X = rW.X , rW.Z = rW.Z )

}



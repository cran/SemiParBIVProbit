working.comp <- function(x, VC = VC, myf = myf){

  e.par <- x$argument
  ll <- length(e.par) 
  p.const <- 2
  
  if(x$PL != "P") {   if(x$eqPL=="both"){ ll <- ll-2; p.const <- 4 } else {ll <- ll-1; p.const <- 3}   } 
  
  e.par <- e.par[-c(ll)]
  good <- x$good
  n <- sum(as.numeric(good==TRUE))

  if(x$PL == "P") X <- rW.X <- Matrix(0,p.const*n,(VC$X1.d2+VC$X2.d2))
  if(x$PL != "P"){ if(x$eqPL == "both") X <- rW.X <- Matrix(0,p.const*n,(VC$X1.d2+VC$X2.d2+2)) else X <- rW.X <- Matrix(0,p.const*n,(VC$X1.d2+VC$X2.d2+1))   } 
  D <- rW.Z <- Matrix(0,p.const*n,1)
  X1 <- VC$X1[good,] 
  X2 <- VC$X2[good,] 

  sqq <- seq(1,(p.const*n-(p.const-1)),by=p.const)

  X[sqq,   1:VC$X1.d2]                      <- X1
  X[sqq+1,(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)] <- X2

  D[sqq,1]   <- x$dl.dbe1
  D[sqq+1,1] <- x$dl.dbe2

  if(x$PL != "P"){ 

    X[sqq+2,VC$X1.d2+VC$X2.d2+1] <- 1 

           if(x$eqPL=="both"){  X[sqq+3,VC$X1.d2+VC$X2.d2+2] <- 1
                                D[sqq+2,1] <- x$dl.dlambda1.st 
                                D[sqq+3,1] <- x$dl.dlambda2.st}  
           if(x$eqPL=="first")  D[sqq+2,1] <- x$dl.dlambda1.st  
           if(x$eqPL=="second") D[sqq+2,1] <- x$dl.dlambda2.st  

                     }


  be1be2 <- cbind(x$d2l.be1.be1, x$d2l.be1.be2, x$d2l.be2.be2)

  if(x$PL != "P"){

       if(x$eqPL=="both")   be1be2 <- cbind(be1be2, x$d2l.be1.lambda1, x$d2l.be1.lambda2, x$d2l.be2.lambda1, x$d2l.be2.lambda2, x$d2l.lambda1.lambda1, 
                                                    x$d2l.lambda1.lambda2, x$d2l.lambda2.lambda2)

       if(x$eqPL=="first")  be1be2 <- cbind(be1be2, x$d2l.be1.lambda1, x$d2l.be2.lambda1, x$d2l.lambda1.lambda1)
       if(x$eqPL=="second") be1be2 <- cbind(be1be2, x$d2l.be1.lambda2, x$d2l.be2.lambda2, x$d2l.lambda2.lambda2)
    
                 }

  L.W <- unlist( apply(be1be2, 1, myf ) , recursive = FALSE) 

  W.inv <- c.W <- list()
  nz <- 0.000000001

    for(i in 1:n){

      W.eig <- eigen(L.W[[i]], symmetric=TRUE)  
      s.ch <- sum(as.numeric(as.vector(L.W[[i]])==0))
      
     if(VC$Model=="BSS" && s.ch==3){ 
      
                          if(L.W[[i]][1,1] < .Machine$double.eps) L.W[[i]][1,1] <- nz
                          c.W[[i]] <- W.inv[[i]] <- matrix(0,2,2)
                          c.W[[i]][1,1] <- sqrt(L.W[[i]][1,1]) 
                          W.inv[[i]][1,1] <- 1/L.W[[i]][1,1]
                          
                                   }else{
                                   
      if(min(W.eig$values) < .Machine$double.eps) { pep <- which(W.eig$values < .Machine$double.eps)
      	                                            W.eig$values[pep] <- nz  
      	                                           }                                   
                                        }
      
      if( !(VC$Model=="BSS" && s.ch==3) ){
      
      c.W[[i]]   <- W.eig$vec%*%tcrossprod(diag(sqrt(W.eig$val)),W.eig$vec) 
      W.inv[[i]] <- W.eig$vec%*%tcrossprod(diag(1/W.eig$val    ),W.eig$vec)
      
      }

                  }


      c.W <- bdiag(c.W)
      W.inv <- bdiag(W.inv)


#WW = nearPD( L.W[[4]], ensureSymmetry = FALSE )$mat
#WW.eig <- eigen(WW, symmetric=TRUE)
#WW.inv <- WW.eig$vec%*%tcrossprod(diag(1/WW.eig$val    ),WW.eig$vec)  
#WW.inv%*%WW
#
#WW.eig <- eigen(L.W[[4]], symmetric=TRUE) 
#WW.eig$values[2] <- 0.0000000001 
#WW = WW.eig$vec%*%tcrossprod(diag(WW.eig$val),WW.eig$vec)
#WW.eig <- eigen(WW, symmetric=TRUE)
#WW.inv <- WW.eig$vec%*%tcrossprod(diag(1/WW.eig$val    ),WW.eig$vec)  
#round(WW.inv%*%WW,5)

## alternative but it is slower ... 
#
#  L.W <- bdiag( unlist( apply(be1be2, 1, myf ) , recursive = FALSE) ) 
#  W.eig <- eigen(L.W, symmetric=TRUE) 
#  c.W   <- Matrix( W.eig$vec%*%tcrossprod(diag(sqrt(W.eig$val)),W.eig$vec) )
#  W.inv <- Matrix( W.eig$vec%*%tcrossprod(diag(1/W.eig$val    ),W.eig$vec) )
#
##################

      rW.X <- as.matrix(c.W%*%X) 
      rW.Z <- as.matrix(c.W%*%( X%*%e.par + W.inv%*%D )) 

rm(c.W, X, W.inv, D)
gc()

                   
 list( rW.X=rW.X , rW.Z=rW.Z )

}



























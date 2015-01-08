working.comp <- function(x, VC = VC, myf = myf,  extra.l =  extra.l){
   
  good <- x$good
  n <- sum(as.numeric(good==TRUE))
  sqq <- seq(1,(extra.l$p.const*n-(extra.l$p.const-1)),by=extra.l$p.const)
  e.par <- x$argument[-c(extra.l$ll)]

  D <- extra.l$D[good,]
  
  D[sqq,1]   <- x$dl.dbe1
  D[sqq+1,1] <- x$dl.dbe2
  
    if(x$PL != "P"){ 
    
             if(x$eqPL=="both"){  D[sqq+2,1] <- x$dl.dlambda1.st 
                                  D[sqq+3,1] <- x$dl.dlambda2.st}  
             if(x$eqPL=="first")  D[sqq+2,1] <- x$dl.dlambda1.st  
             if(x$eqPL=="second") D[sqq+2,1] <- x$dl.dlambda2.st  
    
                 }
  
  
  
  X <- extra.l$X[good,]
  
  X[sqq,   1:VC$X1.d2]                      <- VC$X1[good,] 
  X[sqq+1,(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)] <- VC$X2[good,]   
  
  if(x$PL != "P"){ 
  
      X[sqq+2,VC$X1.d2+VC$X2.d2+1] <- 1 
  
      if(x$eqPL=="both") X[sqq+3,VC$X1.d2+VC$X2.d2+2] <- 1  

                       }


  be1be2 <- cbind(x$d2l.be1.be1, x$d2l.be1.be2, x$d2l.be2.be2)

  if(x$PL != "P"){

       if(x$eqPL=="both")   be1be2 <- cbind(be1be2, x$d2l.be1.lambda1, x$d2l.be1.lambda2, x$d2l.be2.lambda1, x$d2l.be2.lambda2, x$d2l.lambda1.lambda1, 
                                                    x$d2l.lambda1.lambda2, x$d2l.lambda2.lambda2)

       if(x$eqPL=="first")  be1be2 <- cbind(be1be2, x$d2l.be1.lambda1, x$d2l.be2.lambda1, x$d2l.lambda1.lambda1)
       if(x$eqPL=="second") be1be2 <- cbind(be1be2, x$d2l.be1.lambda2, x$d2l.be2.lambda2, x$d2l.lambda2.lambda2)
       
    
                 }

  L.W <- unlist( apply(be1be2, 1, myf ) , recursive = FALSE) 

  W.invsr <- c.W <- list()
  nz <- 0.0000001

    for(i in 1:n){

      
      s.ch <- sum(as.numeric(as.vector(L.W[[i]])==0))
      
     if(VC$Model=="BSS" && s.ch==3){ 
      
                          if(L.W[[i]][1,1] < sqrt(.Machine$double.eps)) L.W[[i]][1,1] <- nz
                          c.W[[i]] <- W.invsr[[i]] <- matrix(0,2,2)
                          c.W[[i]][1,1] <- sqrt(L.W[[i]][1,1]) 
                          W.invsr[[i]][1,1] <- 1/c.W[[i]][1,1]
                          
                                   }else{
                                   
      W.eig <- eigen(L.W[[i]], symmetric=TRUE)                               
      if(min(W.eig$values) < sqrt(.Machine$double.eps)) { pep <- which(W.eig$values < sqrt(.Machine$double.eps))
      	                                            W.eig$values[pep] <- nz  
      	                                           }                                   
                                        }
      
      if( !(VC$Model=="BSS" && s.ch==3) ){
      
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



























bprobgHsPO0 <- function(params, respvec, VC, sp = NULL, qu.mag = NULL){

  

  X1 <- X2 <- 1
  epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)
  max.p   <- 0.9999999

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]

  p1 <- pnorm(eta1); d.n1 <- dnorm(eta1)
  p2 <- pnorm(eta2); d.n2 <- dnorm(eta2) 


  p1 <- ifelse(p1 > max.p , max.p , p1) 
  p1 <- ifelse(p1 < epsilon, epsilon, p1) 
  p2 <- ifelse(p2 > max.p , max.p , p2) 
  p2 <- ifelse(p2 < epsilon, epsilon, p2) 

  X1 <- as.matrix(VC$X1)
  X2 <- as.matrix(VC$X2)

  Y  <- respvec$y1
  CY <- respvec$cy
  
  weights <- VC$weights


########################################################################################################  
  
teta.st <- teta <- 0  
  
p11 <-  pmax(BiCDF(p1, p2, VC$nC, teta), epsilon)

########################################################################################################

  cp11 <- pmax(1 - p11, epsilon)
  
  l.par <- weights*( Y*log(p11) + CY*log(cp11) )


c.copula.be1 <- p2                          
c.copula.be2 <- p1

c.copula2.be1    <- c.copula2.be2 <- 0    
c.copula2.be1be2 <- 1


bit1.b1b1 <- c.copula2.be1*(d.n1)^2-c.copula.be1*d.n1*eta1
bit1.b2b2 <- c.copula2.be2*(d.n2)^2-c.copula.be2*d.n2*eta2
bit1.b1b2 <- c.copula2.be1be2 * d.n1 * d.n2


  dl.dbe1 <-  weights*d.n1*( Y*c.copula.be1/p11   - CY*c.copula.be1/cp11   )                                     
  dl.dbe2 <-  weights*d.n2*( Y*c.copula.be2/p11   - CY*c.copula.be2/cp11   )
                                                  

if(VC$hess==TRUE){

  d2l.be1.be1  <- -weights*(Y*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2 -
                           CY*(bit1.b1b1*cp11+(c.copula.be1*d.n1)^2)/cp11^2 )

  d2l.be2.be2  <- -weights*(Y*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2 -
                           CY*(bit1.b2b2*cp11 + (c.copula.be2*d.n2)^2 )/cp11^2 )

  d2l.be1.be2  <- -weights*(Y*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2 -
                           CY*(bit1.b1b2*cp11+(c.copula.be1*d.n1*c.copula.be2*d.n2))/cp11^2)


}  
  
if(VC$hess==FALSE){
                          
  d2l.be1.be1  <- -weights*((bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11 -
                            (bit1.b1b1*cp11+(c.copula.be1*d.n1)^2)/cp11 )

  d2l.be2.be2  <- -weights*((bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11 -
                            (bit1.b2b2*cp11 + (c.copula.be2*d.n2)^2 )/cp11 )

  d2l.be1.be2  <- -weights*((bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11 -
                            (bit1.b1b2*cp11+(c.copula.be1*d.n1*c.copula.be2*d.n2))/cp11)

                         
}

       
if( is.null(VC$X3) ){

  be1.be1 <- crossprod(X1*c(d2l.be1.be1),X1)
  be2.be2 <- crossprod(X2*c(d2l.be2.be2),X2)
  be1.be2 <- crossprod(X1*c(d2l.be1.be2),X2)

  
  H <- rbind( cbind( be1.be1    , be1.be2 ), 
              cbind( t(be1.be2) , be2.be2 )    
            ) 
            
            
    
         
         G   <- -c( colSums( c(dl.dbe1)*X1 ),
                    colSums( c(dl.dbe2)*X2 )
                  )
    
}




res <- -sum(l.par)

if( ( VC$l.sp1==0 && VC$l.sp2==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)


if(VC$extra.regI == "pC" && VC$hess==FALSE) H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2)  
   
rm(X1, X2)  
     
    

         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps,
              p11=p11, cp11=cp11, eta1=eta1, eta2=eta2, good = rep(TRUE,length(eta1)), 
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2,
              BivD=VC$BivD, p1=p1, p2=p2)      

}




     
























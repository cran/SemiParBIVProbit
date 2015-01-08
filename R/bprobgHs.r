bprobgHs <- function(params,
                     sp.xi1, sp.xi2, 
                     PL, eqPL, valPL, fitPL, 
                     respvec, VC, 
                     sp = NULL, qu.mag = NULL, AT = FALSE){

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]

  p1 <- pnorm(eta1); d.n1 <- dnorm(eta1)
  p2 <- pnorm(eta2); d.n2 <- dnorm(eta2) 

  teta.st <- params[(VC$X1.d2+VC$X2.d2+1)]
  epsilon <- sqrt(.Machine$double.eps)

  criteria <- c(0,1)
  no.good <- apply(apply(cbind(p1,p2), c(1,2), `%in%`, criteria), 1, any)
  good <- no.good==FALSE

  p1 <- p1[good]
  p2 <- p2[good]
  d.n1 <- d.n1[good]
  d.n2 <- d.n2[good]
  eta1 <- eta1[good] 
  eta2 <- eta2[good] 
  X1 <- VC$X1[good,]
  X2 <- VC$X2[good,]
  y1 <- respvec$y1[good]
  y2 <- respvec$y2[good]
  y1.y2 <- respvec$y1.y2[good]
  y1.cy2 <- respvec$y1.cy2[good]
  cy1.y2 <- respvec$cy1.y2[good]
  cy1.cy2 <- respvec$cy1.cy2[good]
  weights <- VC$weights[good]
  
  
  
########################################################################################################  
  
  if( VC$BivD %in% c("N","T") ) teta.st <- ifelse(abs(teta.st) > 8.75, sign(teta.st)*8.75, teta.st )  
  
  if( VC$BivD %in%  c("C0","C180","C90","C270","J0","J180","J90","J270","G0","G180","G90","G270") ) {
 
  if( sign(teta.st)==1  ) teta.st <- ifelse( teta.st > 709, 709, teta.st )  
  if( sign(teta.st)==-1 ) teta.st <- ifelse( teta.st < -20, -20, teta.st )  
  
  }
  
 
    if(VC$BivD %in% c("N","T")      ) teta <- tanh(teta.st) # ; if(teta %in% c(-1,1)) teta <- sign(teta)*0.9999999}
    if(VC$BivD=="F")                  teta <- teta.st + epsilon
    if(VC$BivD %in% c("C0", "C180") ) teta <-    exp(teta.st)       # + epsilon
    if(VC$BivD %in% c("C90","C270") ) teta <- -( exp(teta.st) )     # + epsilon ) 
    if(VC$BivD %in% c("J0", "J180","G0", "G180") ) teta <-    exp(teta.st) + 1   #+ epsilon 
    if(VC$BivD %in% c("J90","J270","G90","G270") ) teta <- -( exp(teta.st) + 1 ) #+ epsilon ) 
    
    #if(VC$BivD %in% c("G0", "G180") ) teta <-    exp(teta.st) + 1 
    #if(VC$BivD %in% c("G90","G270") ) teta <- -( exp(teta.st) + 1 ) 

if(VC$BivD=="N") C.copula <- pbinorm( eta1, eta2, cov12=teta) else C.copula <- BiCopCDF(p1,p2, VC$nC, par=teta, par2=VC$nu)
########################################################################################################


  p11 <- pmax( C.copula, epsilon )
  p10 <- pmax( p1 - p11, epsilon )
  p01 <- pmax( p2 - p11, epsilon )
  p00 <- pmax( 1- p11 - p10 - p01, epsilon )

  l.par <- weights*( y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) )

dH <- copgHs(p1,p2,eta1=NULL,eta2=NULL,teta,teta.st,xi1=NULL,xi1.st=NULL,xi2=NULL,xi2.st=NULL,VC$BivD,VC$nC,VC$nu,PL,eqPL)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
c.copula.theta <- dH$c.copula.theta 
  
c.copula2.be1 <- dH$c.copula2.be1   
c.copula2.be2 <- dH$c.copula2.be2 


bit1.b1b1 <- c.copula2.be1*(d.n1)^2-c.copula.be1*d.n1*eta1
bit2.b1b1 <- -d.n1*eta1-bit1.b1b1
bit3.b1b1 <- -bit1.b1b1
bit4.b1b1 <- -bit2.b1b1

bit1.b2b2 <- c.copula2.be2*(d.n2)^2-c.copula.be2*d.n2*eta2
bit2.b2b2 <- -bit1.b2b2
bit3.b2b2 <- -d.n2*eta2-bit1.b2b2
bit4.b2b2 <- -bit3.b2b2

c.copula2.be1be2 <- dH$c.copula2.be1be2
bit1.b1b2 <- c.copula2.be1be2 * d.n1 *d.n2
bit2.b1b2 <- -bit1.b1b2
bit3.b1b2 <- -bit1.b1b2
bit4.b1b2 <- bit1.b1b2

c.copula2.be1th <- dH$c.copula2.be1th 
bit1.b1th <- c.copula2.be1th*d.n1
bit2.b1th <- -bit1.b1th 
bit3.b1th <- -bit1.b1th 
bit4.b1th <- bit1.b1th 

c.copula2.be2th <- dH$c.copula2.be2th
bit1.b2th <- c.copula2.be2th*d.n2
bit2.b2th <- -bit1.b2th 
bit3.b2th <- -bit1.b2th 
bit4.b2th <- bit1.b2th 

if(AT==TRUE) bit1.th2 <- dH$bit1.th2ATE else bit1.th2 <- dH$bit1.th2

bit2.th2 <- -bit1.th2 
bit3.th2 <- -bit1.th2 
bit4.th2 <- bit1.th2 

  dl.dbe1 <-  weights*d.n1*( (y1.y2*c.copula.be1/p11)  +
                      (y1.cy2*(1-c.copula.be1)/p10) +
                      (cy1.y2*c.copula.be1/(-p01)) +
                      (cy1.cy2*(c.copula.be1-1)/p00) )
                                
  dl.dbe2 <-  weights*d.n2*( (y1.y2*c.copula.be2/p11)  +
                            (y1.cy2*c.copula.be2/(-p10)) +
                                (cy1.y2*(1-c.copula.be2)/(p01)) +
                                (cy1.cy2*(c.copula.be2-1)/p00) )

  dl.drho <- weights*( y1.y2*c.copula.theta/p11+y1.cy2*(-c.copula.theta)/p10 + 
                       cy1.y2*(-c.copula.theta)/p01+cy1.cy2*c.copula.theta/p00 ) 


add.b  <- 1
if(AT==TRUE){
    if(VC$BivD %in% c("N","T") )                                add.b <- 1/cosh(teta.st)^2
    if(VC$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) add.b <-  exp(teta.st)     
    if(VC$BivD %in% c("C90","C270","J90","J270","G90","G270") ) add.b <- -exp(teta.st)        
}



if(VC$hess==TRUE){

  d2l.be1.be1  <- -weights*(y1.y2*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2+
                              y1.cy2*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2+
                              cy1.y2*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2+
                              cy1.cy2*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2 )

  d2l.be2.be2  <- -weights*(y1.y2*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              y1.cy2*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2+
                              cy1.y2*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2+
                              cy1.cy2*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2 )

  d2l.be1.be2  <- -weights*(y1.y2*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2+
                              y1.cy2*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2+
                              cy1.y2*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2+
                              cy1.cy2*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2 )

  d2l.be1.rho  <- -weights*(y1.y2*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )/add.b

  d2l.be2.rho  <- -weights*(y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )/add.b
                              
  d2l.rho.rho  <- (-weights*(   y1.y2*(bit1.th2*p11-( c.copula.theta/add.b)^2)/p11^2+
                               y1.cy2*(bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10^2+
                               cy1.y2*(bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01^2+
                              cy1.cy2*(bit4.th2*p00-( c.copula.theta/add.b)^2)/p00^2 ) )                            
}     


if(VC$hess==FALSE && VC$end==0){

  d2l.be1.be1  <- -weights*(  (bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11+
                              (bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10+
                              (bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01+
                              (bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00 )

  d2l.be2.be2  <- -weights*((bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11+
                              (bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10+
                              (bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01+
                              (bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00 )

  d2l.be1.be2  <- -weights*((bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11+
                              (bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10+
                              (bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01+
                              (bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00 )

  d2l.be1.rho  <- -weights*((bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11+
                              (bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10+
                              (bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01+
                              (bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00 )/add.b

  d2l.be2.rho  <- -weights*(  (bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11+
                              (bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10+
                              (bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01+
                              (bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00 )/add.b

  d2l.rho.rho  <- -weights*(  (bit1.th2*p11-( c.copula.theta/add.b)^2)/p11+
                              (bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10+
                              (bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01+
                              (bit4.th2*p00-( c.copula.theta/add.b)^2)/p00 )
}  

if(VC$hess==FALSE && (VC$end==1 || VC$end==2)){

if(VC$end==1){

fi <- p11/p1
se <- p10/p1
th <- p01/(1-p1)
fo <- p00/(1-p1)

resp1 <- y1
resp2 <- y1
resp3 <- 1 - y1 
resp4 <- 1 - y1 


}

if(VC$end==2){

fi <- p11/p2
se <- p10/(1-p2)
th <- p01/p2
fo <- p00/(1-p2)

resp1 <- y2
resp2 <- 1 - y2
resp3 <- y2 
resp4 <- 1 - y2

}

  d2l.be1.be1  <- -weights*(  fi*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2*resp1+
                              se*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2*resp2+
                              th*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2*resp3+
                              fo*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2*resp4 )

  d2l.be2.be2  <- -weights*(  fi*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2*resp1+
                              se*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2*resp2+
                              th*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2*resp3+
                              fo*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2*resp4 )

  d2l.be1.be2  <- -weights*(  fi*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2*resp1+
                              se*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2*resp2+
                              th*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2*resp3+
                              fo*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2*resp4 )

  d2l.be1.rho  <- -weights*(  fi*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2*resp1+
                              se*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2*resp2+
                              th*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2*resp3+
                              fo*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2*resp4 )/add.b

  d2l.be2.rho  <- -weights*(  fi*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2*resp1+
                              se*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2*resp2+
                              th*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2*resp3+
                              fo*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2*resp4 )/add.b

  d2l.rho.rho  <- -weights*(  fi*(bit1.th2*p11-( c.copula.theta/add.b)^2)/p11^2*resp1+
                              se*(bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10^2*resp2+
                              th*(bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01^2*resp3+
                              fo*(bit4.th2*p00-( c.copula.theta/add.b)^2)/p00^2*resp4 )
}






  be1.be1 <- crossprod(X1*c(d2l.be1.be1),X1)
  be2.be2 <- crossprod(X2*c(d2l.be2.be2),X2)
  be1.be2 <- crossprod(X1*c(d2l.be1.be2),X2)
  be1.rho <- t(t(rowSums(t(X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(X2*c(d2l.be2.rho)))))
  
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
            ) 
            
            
         res <- -sum(l.par)
         
         G   <- -c( colSums( c(dl.dbe1)*X1 ),
                    colSums( c(dl.dbe2)*X2 ),
                    sum( dl.drho )  )
    

if( ( VC$l.sp1==0 && VC$l.sp2==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)


  if(VC$extra.regI == "pC" && VC$hess==FALSE){

      op <- options(warn = -1)
      R <- chol(H, pivot = TRUE)
      options(op)
      p <- dim(H)[2]
      ipiv <- piv <- attr(R, "pivot")
      ipiv[piv] <- 1:p
      rank <- attr(R, "rank")
      ind <- 1:rank
      if (rank < p) R[(rank + 1):p, ] <- 0
      R <- R[ipiv, ipiv]
      H <- crossprod(R)
  
  } 

         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
        
        
  if(VC$extra.regI == "sED"){
  
  ds <- 0
  
  eH <- eigen(H, symmetric = TRUE)
  if(min(eH$values) < epsilon){ eH$values <- abs(eH$values); ds <- 1 }
  if(min(eH$values) < epsilon){ eH$values[which(eH$values < epsilon)] <- 0.0000001; ds <- 1 }
  
  if(ds == 1) H <- eH$vectors%*%tcrossprod(diag(1/eH$values),eH$vectors)   
  
  }  
  

         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps, 
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho,good=good, 
              PL=PL, eqPL=eqPL, BivD=VC$BivD, p1=p1, p2=p2)      

}




     
























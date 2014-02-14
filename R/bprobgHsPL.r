bprobgHsPL <- function(params, BivD, nC, nu, sp.xi1, sp.xi2, PL, eqPL, H.n, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, weights=weights, X1.d2, X2.d2, pPen1=NULL, pPen2=NULL, sp=NULL, qu.mag=NULL, gp1, gp2, fp, l.sp1, l.sp2){

dl.dlambda1.st <- dl.dlambda2.st <- d2l.be1.lambda1 <- d2l.be1.lambda2 <- d2l.be2.lambda1 <- d2l.be2.lambda2 <- d2l.rho.lambda1 <- d2l.rho.lambda2 <- d2l.lambda1.lambda1 <- d2l.lambda2.lambda2 <- d2l.lambda1.lambda2 <- NA 

  epsilon <- .Machine$double.eps*10^6

  eta1 <- X1%*%params[1:X1.d2]
  eta2 <- X2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  teta.st    <- params[(X1.d2+X2.d2+1)]
  
if(eqPL=="both"){  
  lambda1.st <- params[(X1.d2+X2.d2+2)]
  lambda2.st <- params[(X1.d2+X2.d2+3)]
  lambda1 <- exp(lambda1.st)+ epsilon  
  lambda2 <- exp(lambda2.st)+ epsilon
}
if(eqPL=="first"){  
  lambda1.st <- params[(X1.d2+X2.d2+2)]
  lambda2.st <- 0 
  lambda1 <- exp(lambda1.st)+ epsilon  
  lambda2 <- exp(lambda2.st)
}
if(eqPL=="second"){  
  lambda1.st <- 0
  lambda2.st <- params[(X1.d2+X2.d2+2)]
  lambda1 <- exp(lambda1.st) 
  lambda2 <- exp(lambda2.st)+ epsilon
}

  #lambda1 <- xi1 
  #lambda2 <- xi2

  
  if(PL=="PP"){

    if(eqPL=="both"){

    p1 <- pnorm(eta1)^lambda1
    p2 <- pnorm(eta2)^lambda2
    d.n1 <- pnorm(eta1)^(lambda1 - 1) * (lambda1 * dnorm(eta1))
    d.n2 <- pnorm(eta2)^(lambda2 - 1) * (lambda2 * dnorm(eta2))
  
    }
  
    if(eqPL=="first"){
  
    p1 <- pnorm(eta1)^lambda1
    p2 <- pnorm(eta2)
    d.n1 <- pnorm(eta1)^(lambda1 - 1) * (lambda1 * dnorm(eta1))
    d.n2 <- dnorm(eta2)
    
    }
    
    if(eqPL=="second"){
  
    p1 <- pnorm(eta1)
    p2 <- pnorm(eta2)^lambda2
    d.n1 <- dnorm(eta1)
    d.n2 <- pnorm(eta2)^(lambda2 - 1) * (lambda2 * dnorm(eta2))
    
    }    
  
 
  }

  if(PL=="RPP"){


   if(eqPL=="both"){
   
   p1 <- 1-pnorm(-eta1)^lambda1
   p2 <- 1-pnorm(-eta2)^lambda2
   d.n1 <- pnorm(-eta1)^(lambda1 - 1) * (lambda1 * dnorm(-eta1))
   d.n2 <- pnorm(-eta2)^(lambda2 - 1) * (lambda2 * dnorm(-eta2))  
  
   }
   
   if(eqPL=="first"){
   
   p1 <- 1-pnorm(-eta1)^lambda1
   p2 <- pnorm(eta2)
   d.n1 <- pnorm(-eta1)^(lambda1 - 1) * (lambda1 * dnorm(-eta1))
   d.n2 <- dnorm(eta2)  
  
   }   
   
   if(eqPL=="second"){
   
   p1 <- pnorm(eta1)
   p2 <- 1-pnorm(-eta2)^lambda2
   d.n1 <- dnorm(eta1)
   d.n2 <- pnorm(-eta2)^(lambda2 - 1) * (lambda2 * dnorm(-eta2))  
  
   }   
   
  
  }
  

  criteria <- c(0,1)
  no.good <- apply(apply(cbind(p1,p2), c(1,2), `%in%`, criteria), 1, any)
  good <- no.good==FALSE

  p1 <- p1[good]
  p2 <- p2[good]
  d.n1 <- d.n1[good]
  d.n2 <- d.n2[good]
  eta1 <- eta1[good] 
  eta2 <- eta2[good] 
  X1 <- X1[good,]
  X2 <- X2[good,]
  y1.y2 <- y1.y2[good]
  y1.cy2 <- y1.cy2[good]
  cy1.y2 <- cy1.y2[good]
  cy1.cy2 <- cy1.cy2[good]
  weights <- weights[good]

########################################################################################################

    if(BivD %in% c("N","T")      ){teta <- tanh(teta.st); if(teta %in% c(-1,1)) teta <- sign(teta)*0.9999999}
    if(BivD=="F")                  teta <- teta.st + epsilon
    if(BivD %in% c("C0", "C180") ) teta <- exp(teta.st) + epsilon
    if(BivD %in% c("C90","C270") ) teta <- -( exp(teta.st) + epsilon ) 
    if(BivD %in% c("J0", "J180") ) teta <- exp(teta.st) + 1 + epsilon 
    if(BivD %in% c("J90","J270") ) teta <- -( exp(teta.st) + 1 + epsilon ) 
    if(BivD %in% c("G0", "G180") ) teta <- exp(teta.st) + 1 
    if(BivD %in% c("G90","G270") ) teta <- -( exp(teta.st) + 1 ) 

if(BivD=="N") C.copula <- pmax( abs(pbinorm( qnorm(p1), qnorm(p2), cov12=teta)), 1000*.Machine$double.eps ) else C.copula <- BiCopCDF(p1,p2, nC, par=teta, par2=nu)
########################################################################################################


  p11 <- pmax( C.copula, 1000*.Machine$double.eps )
  p10 <- pmax( p1 - p11, 1000*.Machine$double.eps )
  p01 <- pmax( p2 - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )


  l.par <- weights*( y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) )


dH <- copgHs(p1,p2,eta1,eta2,teta,teta.st,xi1=lambda1,xi1.st=lambda1.st,xi2=lambda2,xi2.st=lambda2.st,BivD,nC,nu,PL,eqPL)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
c.copula.theta <- dH$c.copula.theta 

if(eqPL=="both"){c.copula.lambda1 <- dH$c.copula.lambda1
                 c.copula.lambda2 <- dH$c.copula.lambda2}
if(eqPL=="first")  c.copula.lambda1 <- dH$c.copula.lambda1
if(eqPL=="second") c.copula.lambda2 <- dH$c.copula.lambda2


c.copula2.be1    <- dH$c.copula2.be1   
c.copula2.be2    <- dH$c.copula2.be2 
c.copula2.be1be2 <- dH$c.copula2.be1be2
c.copula2.be1th  <- dH$c.copula2.be1th 
c.copula2.be2th  <- dH$c.copula2.be2th
bit1.th2         <- dH$bit1.th2
der.d.n1.be1     <- dH$der.d.n1.be1     
der.d.n2.be2     <- dH$der.d.n2.be2  


if(eqPL=="both"){
bit1.lambda1.2       <- dH$bit1.lambda1.2
bit1.lambda2.2       <- dH$bit1.lambda2.2
c.copula2.be1lambda1 <- dH$c.copula2.be1lambda1
c.copula2.be2lambda2 <- dH$c.copula2.be2lambda2
c.copula2.be1lambda2 <- dH$c.copula2.be1lambda2
c.copula2.be2lambda1 <- dH$c.copula2.be2lambda1
bit1.thlambda1       <- dH$bit1.thlambda1
bit1.thlambda2       <- dH$bit1.thlambda2
bit1.lambda1lambda2  <- dH$bit1.lambda1lambda2

der.p1.lambda1   <- dH$der.p1.lambda1 
der.p2.lambda2   <- dH$der.p2.lambda2 
der2.p1.lambda1  <- dH$der2.p1.lambda1  
der2.p2.lambda2  <- dH$der2.p2.lambda2    
der.d.n1.lambda1 <- dH$der.d.n1.lambda1
der.d.n2.lambda2 <- dH$der.d.n2.lambda2
}

if(eqPL=="first"){
bit1.lambda1.2       <- dH$bit1.lambda1.2
c.copula2.be1lambda1 <- dH$c.copula2.be1lambda1
c.copula2.be2lambda1 <- dH$c.copula2.be2lambda1
bit1.thlambda1       <- dH$bit1.thlambda1
der.p1.lambda1   <- dH$der.p1.lambda1 
der2.p1.lambda1  <- dH$der2.p1.lambda1  
der.d.n1.lambda1 <- dH$der.d.n1.lambda1

}

if(eqPL=="second"){
bit1.lambda2.2       <- dH$bit1.lambda2.2
c.copula2.be2lambda2 <- dH$c.copula2.be2lambda2
c.copula2.be1lambda2 <- dH$c.copula2.be1lambda2
bit1.thlambda2       <- dH$bit1.thlambda2
der.p2.lambda2   <- dH$der.p2.lambda2 
der2.p2.lambda2  <- dH$der2.p2.lambda2    
der.d.n2.lambda2 <- dH$der.d.n2.lambda2
}





bit1.b1b1 <- c.copula2.be1*(d.n1)^2+c.copula.be1*der.d.n1.be1 
bit2.b1b1 <-  der.d.n1.be1-bit1.b1b1
bit3.b1b1 <- -bit1.b1b1
bit4.b1b1 <- -bit2.b1b1

bit1.b2b2 <- c.copula2.be2*(d.n2)^2+c.copula.be2*der.d.n2.be2
bit2.b2b2 <- -bit1.b2b2
bit3.b2b2 <- der.d.n2.be2-bit1.b2b2
bit4.b2b2 <- -bit3.b2b2

bit1.b1b2 <- c.copula2.be1be2 * d.n1 *d.n2
bit2.b1b2 <- -bit1.b1b2
bit3.b1b2 <- -bit1.b1b2
bit4.b1b2 <- bit1.b1b2

bit1.b1th <- c.copula2.be1th*d.n1
bit2.b1th <- -bit1.b1th 
bit3.b1th <- -bit1.b1th 
bit4.b1th <- bit1.b1th 

bit1.b2th <- c.copula2.be2th*d.n2
bit2.b2th <- -bit1.b2th 
bit3.b2th <- -bit1.b2th 
bit4.b2th <- bit1.b2th 

bit2.th2 <- -bit1.th2 
bit3.th2 <- -bit1.th2 
bit4.th2 <- bit1.th2 

if(eqPL=="both"){  

bit2.lambda1.2 <- der2.p1.lambda1-bit1.lambda1.2
bit3.lambda1.2 <- -bit1.lambda1.2
bit4.lambda1.2<- -bit2.lambda1.2
        
bit2.lambda2.2 <- -bit1.lambda2.2
bit3.lambda2.2 <- der2.p2.lambda2-bit1.lambda2.2
bit4.lambda2.2<- -bit3.lambda2.2
    
bit1.be1lambda1 <- c.copula2.be1lambda1*d.n1+c.copula.be1*der.d.n1.lambda1
bit2.be1lambda1 <- der.d.n1.lambda1- bit1.be1lambda1
bit3.be1lambda1 <-  -bit1.be1lambda1
bit4.be1lambda1 <-  -bit2.be1lambda1

bit1.be2lambda2 <- c.copula2.be2lambda2*d.n2+c.copula.be2*der.d.n2.lambda2
bit2.be2lambda2 <- -bit1.be2lambda2 
bit3.be2lambda2 <- der.d.n2.lambda2 - bit1.be2lambda2
bit4.be2lambda2 <- -bit3.be2lambda2 

bit1.be1lambda2 <- c.copula2.be1lambda2*d.n1
bit2.be1lambda2 <- -bit1.be1lambda2
bit3.be1lambda2 <- -bit1.be1lambda2
bit4.be1lambda2 <- bit1.be1lambda2

bit1.be2lambda1 <- c.copula2.be2lambda1*d.n2
bit2.be2lambda1 <- -bit1.be2lambda1
bit3.be2lambda1 <- -bit1.be2lambda1
bit4.be2lambda1 <- bit1.be2lambda1

bit2.thlambda1 <- -bit1.thlambda1
bit3.thlambda1 <- -bit1.thlambda1
bit4.thlambda1 <- bit1.thlambda1

bit2.thlambda2 <- -bit1.thlambda2
bit3.thlambda2 <- -bit1.thlambda2
bit4.thlambda2 <- bit1.thlambda2

bit2.lambda1lambda2 <- -bit1.lambda1lambda2 
bit3.lambda1lambda2 <- -bit1.lambda1lambda2 
bit4.lambda1lambda2 <- bit1.lambda1lambda2 
}

if(eqPL=="first"){     
bit2.lambda1.2 <- der2.p1.lambda1-bit1.lambda1.2
bit3.lambda1.2 <- -bit1.lambda1.2
bit4.lambda1.2<- -bit2.lambda1.2
       
bit1.be1lambda1 <- c.copula2.be1lambda1*d.n1+c.copula.be1*der.d.n1.lambda1
bit2.be1lambda1 <- der.d.n1.lambda1- bit1.be1lambda1
bit3.be1lambda1 <-  -bit1.be1lambda1
bit4.be1lambda1 <-  -bit2.be1lambda1

bit1.be2lambda1 <- c.copula2.be2lambda1*d.n2
bit2.be2lambda1 <- -bit1.be2lambda1
bit3.be2lambda1 <- -bit1.be2lambda1
bit4.be2lambda1 <- bit1.be2lambda1

bit2.thlambda1 <- -bit1.thlambda1
bit3.thlambda1 <- -bit1.thlambda1
bit4.thlambda1 <- bit1.thlambda1

}

if(eqPL=="second"){     
  
bit2.lambda2.2 <- -bit1.lambda2.2
bit3.lambda2.2 <- der2.p2.lambda2-bit1.lambda2.2
bit4.lambda2.2<- -bit3.lambda2.2

bit1.be2lambda2 <- c.copula2.be2lambda2*d.n2+c.copula.be2*der.d.n2.lambda2
bit2.be2lambda2 <- -bit1.be2lambda2 
bit3.be2lambda2 <- der.d.n2.lambda2 - bit1.be2lambda2
bit4.be2lambda2 <- -bit3.be2lambda2 

bit1.be1lambda2 <- c.copula2.be1lambda2*d.n1
bit2.be1lambda2 <- -bit1.be1lambda2
bit3.be1lambda2 <- -bit1.be1lambda2
bit4.be1lambda2 <- bit1.be1lambda2

bit2.thlambda2 <- -bit1.thlambda2
bit3.thlambda2 <- -bit1.thlambda2
bit4.thlambda2 <- bit1.thlambda2

}



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

if(eqPL=="both"){ 

  dl.dlambda1.st <- weights*(y1.y2*c.copula.lambda1/p11+y1.cy2*(der.p1.lambda1-c.copula.lambda1)/p10 + 
                       cy1.y2*(-c.copula.lambda1)/p01+cy1.cy2*(c.copula.lambda1-der.p1.lambda1)/p00)  
  
  dl.dlambda2.st <- weights*(y1.y2*c.copula.lambda2/p11+y1.cy2*(-c.copula.lambda2)/p10 + 
                       cy1.y2*(der.p2.lambda2-c.copula.lambda2)/p01+cy1.cy2*(c.copula.lambda2-der.p2.lambda2)/p00)

}

if(eqPL=="first"){ 

  dl.dlambda1.st <- weights*(y1.y2*c.copula.lambda1/p11+y1.cy2*(der.p1.lambda1-c.copula.lambda1)/p10 + 
                       cy1.y2*(-c.copula.lambda1)/p01+cy1.cy2*(c.copula.lambda1-der.p1.lambda1)/p00)  

}

if(eqPL=="second"){  
  
  dl.dlambda2.st <- weights*(y1.y2*c.copula.lambda2/p11+y1.cy2*(-c.copula.lambda2)/p10 + 
                       cy1.y2*(der.p2.lambda2-c.copula.lambda2)/p01+cy1.cy2*(c.copula.lambda2-der.p2.lambda2)/p00)

}


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
                              cy1.cy2*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )

  d2l.be2.rho  <- -weights*(y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )

  d2l.rho.rho  <- -weights*(y1.y2*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                              y1.cy2*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2+
                              cy1.y2*(bit3.th2*p01-(-c.copula.theta)^2)/p01^2+
                              cy1.cy2*(bit4.th2*p00-c.copula.theta^2)/p00^2 )


if(eqPL=="both"){  
                           
  d2l.be1.lambda1  <- -weights*(y1.y2*(bit1.be1lambda1*p11-(c.copula.be1*d.n1*c.copula.lambda1))/p11^2+
                              y1.cy2*(bit2.be1lambda1*p10-((1-c.copula.be1)*d.n1)*(der.p1.lambda1-c.copula.lambda1))/p10^2+
                              cy1.y2*(bit3.be1lambda1*p01-(-c.copula.be1*d.n1)*(-c.copula.lambda1))/p01^2+
                              cy1.cy2*(bit4.be1lambda1*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda1-der.p1.lambda1))/p00^2 )  
  
  d2l.be1.lambda2  <-  -weights*(y1.y2*(bit1.be1lambda2*p11-(c.copula.be1*d.n1*c.copula.lambda2))/p11^2+
                              y1.cy2*(bit2.be1lambda2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.lambda2))/p10^2+
                              cy1.y2*(bit3.be1lambda2*p01-(-c.copula.be1*d.n1)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              cy1.cy2*(bit4.be1lambda2*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda2-der.p2.lambda2))/p00^2  )
  
  d2l.be2.lambda1  <-  -weights*(y1.y2*(bit1.be2lambda1*p11-(c.copula.be2*d.n2*c.copula.lambda1))/p11^2+
                              y1.cy2*(bit2.be2lambda1*p10-(-c.copula.be2*d.n2)*(der.p1.lambda1-c.copula.lambda1))/p10^2+
                              cy1.y2*(bit3.be2lambda1*p01-((1-c.copula.be2)*d.n2)*(-c.copula.lambda1))/p01^2+
                              cy1.cy2*(bit4.be2lambda1*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda1-der.p1.lambda1))/p00^2   ) 
    
  d2l.be2.lambda2  <- -weights*(y1.y2*(bit1.be2lambda2*p11-(c.copula.be2*d.n2*c.copula.lambda2))/p11^2+
                              y1.cy2*(bit2.be2lambda2*p10-(-c.copula.be2*d.n2)*(-c.copula.lambda2))/p10^2+
                              cy1.y2*(bit3.be2lambda2*p01-((1-c.copula.be2)*d.n2)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              cy1.cy2*(bit4.be2lambda2*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda2-der.p2.lambda2))/p00^2 )  
  
  d2l.rho.lambda1  <- -weights*(y1.y2*(bit1.thlambda1*p11-c.copula.theta*c.copula.lambda1)/p11^2+
                              y1.cy2*(bit2.thlambda1*p10-(-c.copula.theta*(der.p1.lambda1-c.copula.lambda1)))/p10^2+
                              cy1.y2*(bit3.thlambda1*p01-(-c.copula.theta*(-c.copula.lambda1)))/p01^2+
                              cy1.cy2*(bit4.thlambda1*p00-c.copula.theta*(c.copula.lambda1-der.p1.lambda1))/p00^2  )
    
  d2l.rho.lambda2  <- -weights*(y1.y2*(bit1.thlambda2*p11-c.copula.theta*c.copula.lambda2)/p11^2+
                              y1.cy2*(bit2.thlambda2*p10-(-c.copula.theta*(-c.copula.lambda2)))/p10^2+
                              cy1.y2*(bit3.thlambda2*p01-(-c.copula.theta*(der.p2.lambda2-c.copula.lambda2)))/p01^2+
                              cy1.cy2*(bit4.thlambda2*p00-c.copula.theta*(c.copula.lambda2-der.p2.lambda2))/p00^2  )
  
  d2l.lambda1.lambda1  <-   -weights*(y1.y2*(bit1.lambda1.2*p11-c.copula.lambda1^2)/p11^2+
                              y1.cy2*(bit2.lambda1.2*p10-(der.p1.lambda1-c.copula.lambda1)^2)/p10^2+
                              cy1.y2*(bit3.lambda1.2*p01-(-c.copula.lambda1)^2)/p01^2+
                              cy1.cy2*(bit4.lambda1.2*p00-(c.copula.lambda1-der.p1.lambda1)^2)/p00^2)
    
  d2l.lambda2.lambda2  <- -weights*(y1.y2*(bit1.lambda2.2*p11-c.copula.lambda2^2)/p11^2+
                              y1.cy2*(bit2.lambda2.2*p10-(-c.copula.lambda2)^2)/p10^2+
                              cy1.y2*(bit3.lambda2.2*p01-(der.p2.lambda2-c.copula.lambda2)^2)/p01^2+
                              cy1.cy2*(bit4.lambda2.2*p00-(c.copula.lambda2-der.p2.lambda2)^2)/p00^2)
  
  d2l.lambda1.lambda2  <- -weights*(y1.y2*(bit1.lambda1lambda2*p11-c.copula.lambda1*c.copula.lambda2)/p11^2+
                              y1.cy2*(bit2.lambda1lambda2*p10-(der.p1.lambda1-c.copula.lambda1)*(-c.copula.lambda2))/p10^2+
                              cy1.y2*(bit3.lambda1lambda2*p01-(-c.copula.lambda1)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              cy1.cy2*(bit4.lambda1lambda2*p00-(c.copula.lambda1-der.p1.lambda1)*(c.copula.lambda2-der.p2.lambda2))/p00^2)


}

if(eqPL=="first"){  
                           
  d2l.be1.lambda1  <- -weights*(y1.y2*(bit1.be1lambda1*p11-(c.copula.be1*d.n1*c.copula.lambda1))/p11^2+
                              y1.cy2*(bit2.be1lambda1*p10-((1-c.copula.be1)*d.n1)*(der.p1.lambda1-c.copula.lambda1))/p10^2+
                              cy1.y2*(bit3.be1lambda1*p01-(-c.copula.be1*d.n1)*(-c.copula.lambda1))/p01^2+
                              cy1.cy2*(bit4.be1lambda1*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda1-der.p1.lambda1))/p00^2 )  
  
  d2l.be2.lambda1  <-  -weights*(y1.y2*(bit1.be2lambda1*p11-(c.copula.be2*d.n2*c.copula.lambda1))/p11^2+
                              y1.cy2*(bit2.be2lambda1*p10-(-c.copula.be2*d.n2)*(der.p1.lambda1-c.copula.lambda1))/p10^2+
                              cy1.y2*(bit3.be2lambda1*p01-((1-c.copula.be2)*d.n2)*(-c.copula.lambda1))/p01^2+
                              cy1.cy2*(bit4.be2lambda1*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda1-der.p1.lambda1))/p00^2   ) 
    
  d2l.rho.lambda1  <- -weights*(y1.y2*(bit1.thlambda1*p11-c.copula.theta*c.copula.lambda1)/p11^2+
                              y1.cy2*(bit2.thlambda1*p10-(-c.copula.theta*(der.p1.lambda1-c.copula.lambda1)))/p10^2+
                              cy1.y2*(bit3.thlambda1*p01-(-c.copula.theta*(-c.copula.lambda1)))/p01^2+
                              cy1.cy2*(bit4.thlambda1*p00-c.copula.theta*(c.copula.lambda1-der.p1.lambda1))/p00^2  )
    
  d2l.lambda1.lambda1  <-   -weights*(y1.y2*(bit1.lambda1.2*p11-c.copula.lambda1^2)/p11^2+
                              y1.cy2*(bit2.lambda1.2*p10-(der.p1.lambda1-c.copula.lambda1)^2)/p10^2+
                              cy1.y2*(bit3.lambda1.2*p01-(-c.copula.lambda1)^2)/p01^2+
                              cy1.cy2*(bit4.lambda1.2*p00-(c.copula.lambda1-der.p1.lambda1)^2)/p00^2)
    

}


if(eqPL=="second"){                             
 
  
  d2l.be1.lambda2  <-  -weights*(y1.y2*(bit1.be1lambda2*p11-(c.copula.be1*d.n1*c.copula.lambda2))/p11^2+
                              y1.cy2*(bit2.be1lambda2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.lambda2))/p10^2+
                              cy1.y2*(bit3.be1lambda2*p01-(-c.copula.be1*d.n1)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              cy1.cy2*(bit4.be1lambda2*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda2-der.p2.lambda2))/p00^2  )


    
  d2l.be2.lambda2  <- -weights*(y1.y2*(bit1.be2lambda2*p11-(c.copula.be2*d.n2*c.copula.lambda2))/p11^2+
                              y1.cy2*(bit2.be2lambda2*p10-(-c.copula.be2*d.n2)*(-c.copula.lambda2))/p10^2+
                              cy1.y2*(bit3.be2lambda2*p01-((1-c.copula.be2)*d.n2)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              cy1.cy2*(bit4.be2lambda2*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda2-der.p2.lambda2))/p00^2 )  


    
  d2l.rho.lambda2  <- -weights*(y1.y2*(bit1.thlambda2*p11-c.copula.theta*c.copula.lambda2)/p11^2+
                              y1.cy2*(bit2.thlambda2*p10-(-c.copula.theta*(-c.copula.lambda2)))/p10^2+
                              cy1.y2*(bit3.thlambda2*p01-(-c.copula.theta*(der.p2.lambda2-c.copula.lambda2)))/p01^2+
                              cy1.cy2*(bit4.thlambda2*p00-c.copula.theta*(c.copula.lambda2-der.p2.lambda2))/p00^2  )


    
  d2l.lambda2.lambda2  <- -weights*(y1.y2*(bit1.lambda2.2*p11-c.copula.lambda2^2)/p11^2+
                              y1.cy2*(bit2.lambda2.2*p10-(-c.copula.lambda2)^2)/p10^2+
                              cy1.y2*(bit3.lambda2.2*p01-(der.p2.lambda2-c.copula.lambda2)^2)/p01^2+
                              cy1.cy2*(bit4.lambda2.2*p00-(c.copula.lambda2-der.p2.lambda2)^2)/p00^2)
  

}


  be1.be1 <- crossprod(X1*c(d2l.be1.be1),X1)
  be2.be2 <- crossprod(X2*c(d2l.be2.be2),X2)
  be1.be2 <- crossprod(X1*c(d2l.be1.be2),X2)
  be1.rho <- t(t(rowSums(t(X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(X2*c(d2l.be2.rho)))))
  
if(eqPL=="both"){ 
  be1.lambda1 <- t(t(rowSums(t(X1*c(d2l.be1.lambda1)))))
  be2.lambda2 <- t(t(rowSums(t(X2*c(d2l.be2.lambda2)))))
  be1.lambda2 <- t(t(rowSums(t(X1*c(d2l.be1.lambda2)))))
  be2.lambda1 <- t(t(rowSums(t(X2*c(d2l.be2.lambda1)))))
}
if(eqPL=="first"){ 
  be1.lambda1 <- t(t(rowSums(t(X1*c(d2l.be1.lambda1)))))
  be2.lambda1 <- t(t(rowSums(t(X2*c(d2l.be2.lambda1)))))
}
if(eqPL=="second"){ 
  be2.lambda2 <- t(t(rowSums(t(X2*c(d2l.be2.lambda2)))))
  be1.lambda2 <- t(t(rowSums(t(X1*c(d2l.be1.lambda2)))))
}




  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) )
            ) 

         G   <- c( colSums( c(dl.dbe1)*X1 ),
                   colSums( c(dl.dbe2)*X2 ),
                   sum( dl.drho )  ) 





if(eqPL=="both"){  


lbit1 <- rbind( cbind( t(be1.lambda1) , t(be2.lambda1) , sum(d2l.rho.lambda1) ),
                cbind( t(be1.lambda2) , t(be2.lambda2) , sum(d2l.rho.lambda2) )
              )
             
lbit2 <- rbind( cbind( be1.lambda1, be1.lambda2 ),
                cbind( be2.lambda1, be2.lambda2 ),
                cbind( sum(d2l.rho.lambda1), sum(d2l.rho.lambda2) ),
                cbind( sum(d2l.lambda1.lambda1), sum(d2l.lambda1.lambda2) ),
                cbind( sum(d2l.lambda1.lambda2), sum(d2l.lambda2.lambda2) )
             )  
             
H <- cbind( rbind( H, lbit1), lbit2 )
G <- -c(G, sum( dl.dlambda1.st ), sum( dl.dlambda2.st ) )

}

if(eqPL=="first"){  


lbit1 <- cbind( t(be1.lambda1) , t(be2.lambda1) , sum(d2l.rho.lambda1) )

lbit2 <- rbind( be1.lambda1,
                be2.lambda1,
                sum(d2l.rho.lambda1),
                sum(d2l.lambda1.lambda1) 
             )  
             
H <- cbind( rbind( H, lbit1), lbit2 )
G <- -c(G, sum( dl.dlambda1.st ) )

}


if(eqPL=="second"){  

lbit1 <- cbind( t(be1.lambda2) , t(be2.lambda2) , sum(d2l.rho.lambda2) )

lbit2 <- rbind( be1.lambda2,
                be2.lambda2,
                sum(d2l.rho.lambda2),
                sum(d2l.lambda2.lambda2) 
             )  
             
H <- cbind( rbind( H, lbit1), lbit2 )
G <- -c(G, sum( dl.dlambda2.st ) )

}

    res <- -sum(l.par)

    if(eqPL=="both"){   add.z <- diag(c(1,1)); add.z <- add.z*c(sp.xi1,sp.xi2); nsh <- 2  } 
    if(eqPL=="first"){  add.z <- sp.xi1; nsh <- 1}
    if(eqPL=="second"){ add.z <- sp.xi2; nsh <- 1}


if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE){ 
    lps <- length(params) - nsh 
    S.h <- adiag( matrix(0,lps,lps), add.z) 
                                         }else{
        
    dimP1 <- dimP2 <- 0     
    S1 <- S2 <- matrix(0,1,1)  

    S <- mapply("*", qu.mag$Ss[-qu.mag$exclu], sp, SIMPLIFY=FALSE)
    S <- do.call(adiag, lapply(S, unlist))

    ma1 <- matrix(0,gp1,gp1) 
    ma2 <- matrix(0,gp2,gp2)

    if(length(pPen1)!=0){ indP1 <- qu.mag$off[1]:(qu.mag$off[1]+qu.mag$rank[1]-1)
                          dimP1 <- length(indP1)
                          ma1[indP1,indP1] <- S[1:dimP1,1:dimP1]
                                } 

    if(length(pPen2)!=0){ 
                          indP2 <- (qu.mag$off[l.sp1+1]-X1.d2):(-X1.d2+qu.mag$off[l.sp1+1]+qu.mag$rank[l.sp1+1]-1)
                          dimP2 <- length(indP2)
                          ma2[indP2,indP2] <- S[(dimP1+1):(length(indP2)+dimP1),(dimP1+1):(length(indP2)+dimP1)]
                                }                                 
    
    lP1 <- length(pPen1); lP2 <- length(pPen2) 
    
    if((lP1!=0 && l.sp1>1) || (lP1==0 && l.sp1>0)) S1 <- S[(dimP1+1):(dimP1+X1.d2-gp1),(dimP1+1):(dimP1+X1.d2-gp1)]
    if((lP2!=0 && l.sp2>1) || (lP2==0 && l.sp2>0)){dS1 <- dim(S1)[2]; if(dS1==1) dS1 <- 0; 
                                                   S2 <- S[(dimP1+dimP2+dS1+1):dim(S)[2],(dimP1+dimP2+dS1+1):dim(S)[2]]}
    
    lS1 <- length(S1); lS2 <- length(S2) 

    if(lS1==1 && lS2==1) S.h <- adiag(ma1, ma2, 0, add.z)
    if(lS1 >1 && lS2==1) S.h <- adiag(ma1, S1, ma2, 0, add.z)
    if(lS1==1 && lS2 >1) S.h <- adiag(ma1, ma2, S2, 0, add.z)
    if(lS1 >1 && lS2 >1) S.h <- adiag(ma1, S1, ma2, S2, 0, add.z)
        
         }

   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params


         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
         H   <- H + S.h  

         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, l.par=l.par, 
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho,
              dl.dlambda1.st=dl.dlambda1.st, dl.dlambda2.st=dl.dlambda2.st,
              d2l.be1.lambda1=d2l.be1.lambda1, d2l.be1.lambda2=d2l.be1.lambda2,    
              d2l.be2.lambda1=d2l.be2.lambda1, d2l.be2.lambda2=d2l.be2.lambda2, 
              d2l.rho.lambda1=d2l.rho.lambda1, d2l.rho.lambda2=d2l.rho.lambda2, 
              d2l.lambda1.lambda1=d2l.lambda1.lambda1, d2l.lambda2.lambda2=d2l.lambda2.lambda2,  
              d2l.lambda1.lambda2=d2l.lambda1.lambda2, good=good, PL=PL, eqPL=eqPL,BivD=BivD)      

}




     
























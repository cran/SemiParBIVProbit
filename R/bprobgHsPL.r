bprobgHsPL <- function(params, sp.xi1, sp.xi2, 
                       PL, eqPL, valPL, fitPL, 
                       respvec, VC, 
                       sp = NULL, qu.mag = NULL){

dl.dlambda1.st <- dl.dlambda2.st <- d2l.be1.lambda1 <- d2l.be1.lambda2 <- d2l.be2.lambda1 <- d2l.be2.lambda2 <- d2l.rho.lambda1 <- d2l.rho.lambda2 <- d2l.lambda1.lambda1 <- d2l.lambda2.lambda2 <- d2l.lambda1.lambda2 <- NA 

  epsilon <- .Machine$double.eps*10^6

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  teta.st <- params[(VC$X1.d2+VC$X2.d2+1)]
  
  
if(PL=="PP" || PL=="RPP"){  
  
if(fitPL=="fixed"){

  lambda1.st <- valPL[1]
  lambda2.st <- valPL[2]
  lambda1 <- exp(lambda1.st)+ epsilon  
  lambda2 <- exp(lambda2.st)+ epsilon

}else{
  
if(eqPL=="both"){  
  lambda1.st <- params[(VC$X1.d2+VC$X2.d2+2)]
  lambda2.st <- params[(VC$X1.d2+VC$X2.d2+3)]
  lambda1 <- exp(lambda1.st)+ epsilon  
  lambda2 <- exp(lambda2.st)+ epsilon
}
if(eqPL=="first"){  
  lambda1.st <- params[(VC$X1.d2+VC$X2.d2+2)]
  lambda2.st <- 0 
  lambda1 <- exp(lambda1.st)+ epsilon  
  lambda2 <- exp(lambda2.st)
}
if(eqPL=="second"){  
  lambda1.st <- 0
  lambda2.st <- params[(VC$X1.d2+VC$X2.d2+2)]
  lambda1 <- exp(lambda1.st) 
  lambda2 <- exp(lambda2.st)+ epsilon
}

}


}

if(PL=="SN"){  

if(fitPL=="fixed"){

  lambda1.st <- valPL[1]
  lambda2.st <- valPL[2]
  lambda1 <- lambda1.st
  lambda2 <- lambda2.st
  del1 <- -lambda1/sqrt(1+lambda1^2)
  del2 <- -lambda2/sqrt(1+lambda2^2)

}else{

if(eqPL=="both"){  
  lambda1.st <- params[(VC$X1.d2+VC$X2.d2+2)]
  lambda2.st <- params[(VC$X1.d2+VC$X2.d2+3)]
  lambda1 <- lambda1.st
  lambda2 <- lambda2.st
  del1 <- -lambda1/sqrt(1+lambda1^2)
  del2 <- -lambda2/sqrt(1+lambda2^2)
}
if(eqPL=="first"){  
  lambda1.st <- params[(VC$X1.d2+VC$X2.d2+2)]
  lambda2.st <- 0 
  lambda1 <- lambda1.st
  lambda2 <- lambda2.st
  del1 <- -lambda1/sqrt(1+lambda1^2)

}
if(eqPL=="second"){  
  lambda1.st <- 0
  lambda2.st <- params[(VC$X1.d2+VC$X2.d2+2)]
  lambda1 <- lambda1.st 
  lambda2 <- lambda2.st
  del2 <- -lambda2/sqrt(1+lambda2^2)
}

}


}

  
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
  
  
  
  if(PL=="SN"){  
  
      if(eqPL=="both"){
  
      p1 <- 2*pmax(pbinorm( eta1, 0, cov12=del1), 1000*.Machine$double.eps )
      p2 <- 2*pmax(pbinorm( eta2, 0, cov12=del2), 1000*.Machine$double.eps )
      d.n1 <- 2*dnorm(eta1)*pnorm(lambda1*eta1)
      d.n2 <- 2*dnorm(eta2)*pnorm(lambda2*eta2)
    
      }
    
      if(eqPL=="first"){
    
      p1 <- 2*pmax(pbinorm( eta1, 0, cov12=del1), 1000*.Machine$double.eps )
      p2 <- pnorm(eta2)
      d.n1 <- 2*dnorm(eta1)*pnorm(lambda1*eta1)
      d.n2 <- dnorm(eta2)
      
      }
      
      if(eqPL=="second"){
    
      p1 <- pnorm(eta1)
      p2 <- 2*pmax(pbinorm( eta2, 0, cov12=del2), 1000*.Machine$double.eps )
      d.n1 <- dnorm(eta1)
      d.n2 <- 2*dnorm(eta2)*pnorm(lambda2*eta2)
      
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

    if(VC$BivD %in% c("N","T")      ){teta <- tanh(teta.st); if(teta %in% c(-1,1)) teta <- sign(teta)*0.9999999}
    if(VC$BivD=="F")                  teta <- teta.st + epsilon
    if(VC$BivD %in% c("C0", "C180") ) teta <- exp(teta.st) + epsilon
    if(VC$BivD %in% c("C90","C270") ) teta <- -( exp(teta.st) + epsilon ) 
    if(VC$BivD %in% c("J0", "J180") ) teta <- exp(teta.st) + 1 + epsilon 
    if(VC$BivD %in% c("J90","J270") ) teta <- -( exp(teta.st) + 1 + epsilon ) 
    if(VC$BivD %in% c("G0", "G180") ) teta <- exp(teta.st) + 1 
    if(VC$BivD %in% c("G90","G270") ) teta <- -( exp(teta.st) + 1 ) 

if(VC$BivD=="N") C.copula <- pbinorm( qnorm(p1), qnorm(p2), cov12=teta) else C.copula <- BiCopCDF(p1,p2, VC$nC, par=teta, par2=VC$nu)
########################################################################################################


  p11 <- pmax( C.copula, 1000*.Machine$double.eps )
  p10 <- pmax( p1 - p11, 1000*.Machine$double.eps )
  p01 <- pmax( p2 - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )


  l.par <- weights*( y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) )


dH <- copgHs(p1,p2,eta1,eta2,teta,teta.st,xi1=lambda1,xi1.st=lambda1.st,xi2=lambda2,xi2.st=lambda2.st,VC$BivD,VC$nC,VC$nu,PL,eqPL)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
c.copula.theta <- dH$c.copula.theta 

if(fitPL!="fixed"){

if(eqPL=="both"){c.copula.lambda1 <- dH$c.copula.lambda1
                 c.copula.lambda2 <- dH$c.copula.lambda2}
if(eqPL=="first")  c.copula.lambda1 <- dH$c.copula.lambda1
if(eqPL=="second") c.copula.lambda2 <- dH$c.copula.lambda2

}


c.copula2.be1    <- dH$c.copula2.be1   
c.copula2.be2    <- dH$c.copula2.be2 
c.copula2.be1be2 <- dH$c.copula2.be1be2
c.copula2.be1th  <- dH$c.copula2.be1th 
c.copula2.be2th  <- dH$c.copula2.be2th
bit1.th2         <- dH$bit1.th2
der.d.n1.be1     <- dH$der.d.n1.be1     
der.d.n2.be2     <- dH$der.d.n2.be2  


if(fitPL!="fixed"){

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



if(fitPL!="fixed"){


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

}

################
# GRADIENT
################


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


if(fitPL!="fixed"){



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

}



#################

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
                              cy1.cy2*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )

  d2l.be2.rho  <- -weights*(y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )

  d2l.rho.rho  <- -weights*(y1.y2*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                              y1.cy2*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2+
                              cy1.y2*(bit3.th2*p01-(-c.copula.theta)^2)/p01^2+
                              cy1.cy2*(bit4.th2*p00-c.copula.theta^2)/p00^2 )
                              
}      

if(VC$hess==FALSE && VC$end==0){


  d2l.be1.be1  <- -weights*((bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11+
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
                              (bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00 )

  d2l.be2.rho  <- -weights*((bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11+
                              (bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10+
                              (bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01+
                              (bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00 )

  d2l.rho.rho  <- -weights*((bit1.th2*p11-c.copula.theta^2)/p11+
                              (bit2.th2*p10-(-c.copula.theta)^2)/p10+
                              (bit3.th2*p01-(-c.copula.theta)^2)/p01+
                              (bit4.th2*p00-c.copula.theta^2)/p00 )
                              
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

  d2l.be1.be1  <- -weights*(  resp1*fi*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2+
                              resp2*se*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2+
                              resp3*th*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2+
                              resp4*fo*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2 )

  d2l.be2.be2  <- -weights*(  resp1*fi*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              resp2*se*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2+
                              resp3*th*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2+
                              resp4*fo*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2 )

  d2l.be1.be2  <- -weights*(  resp1*fi*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2+
                              resp2*se*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2+
                              resp3*th*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2+
                              resp4*fo*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2 )

  d2l.be1.rho  <- -weights*(  resp1*fi*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2+
                              resp2*se*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2+
                              resp3*th*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2+
                              resp4*fo*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )

  d2l.be2.rho  <- -weights*(  resp1*fi*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              resp2*se*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                              resp3*th*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              resp4*fo*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )

  d2l.rho.rho  <- -weights*(  resp1*fi*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                              resp2*se*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2+
                              resp3*th*(bit3.th2*p01-(-c.copula.theta)^2)/p01^2+
                              resp4*fo*(bit4.th2*p00-c.copula.theta^2)/p00^2 )
} 














if(VC$hess==TRUE){


if(fitPL!="fixed"){



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


}



}
#################

if(VC$hess==FALSE && VC$end==0){


if(fitPL!="fixed"){



if(eqPL=="both"){  
                           
  d2l.be1.lambda1  <- -weights*((bit1.be1lambda1*p11-(c.copula.be1*d.n1*c.copula.lambda1))/p11+
                              (bit2.be1lambda1*p10-((1-c.copula.be1)*d.n1)*(der.p1.lambda1-c.copula.lambda1))/p10+
                              (bit3.be1lambda1*p01-(-c.copula.be1*d.n1)*(-c.copula.lambda1))/p01+
                              (bit4.be1lambda1*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda1-der.p1.lambda1))/p00 )  
  
  d2l.be1.lambda2  <-  -weights*((bit1.be1lambda2*p11-(c.copula.be1*d.n1*c.copula.lambda2))/p11+
                              (bit2.be1lambda2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.lambda2))/p10+
                              (bit3.be1lambda2*p01-(-c.copula.be1*d.n1)*(der.p2.lambda2-c.copula.lambda2))/p01+
                              (bit4.be1lambda2*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda2-der.p2.lambda2))/p00  )
  
  d2l.be2.lambda1  <-  -weights*((bit1.be2lambda1*p11-(c.copula.be2*d.n2*c.copula.lambda1))/p11+
                              (bit2.be2lambda1*p10-(-c.copula.be2*d.n2)*(der.p1.lambda1-c.copula.lambda1))/p10+
                              (bit3.be2lambda1*p01-((1-c.copula.be2)*d.n2)*(-c.copula.lambda1))/p01+
                              (bit4.be2lambda1*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda1-der.p1.lambda1))/p00   ) 
    
  d2l.be2.lambda2  <- -weights*((bit1.be2lambda2*p11-(c.copula.be2*d.n2*c.copula.lambda2))/p11+
                              (bit2.be2lambda2*p10-(-c.copula.be2*d.n2)*(-c.copula.lambda2))/p10+
                              (bit3.be2lambda2*p01-((1-c.copula.be2)*d.n2)*(der.p2.lambda2-c.copula.lambda2))/p01+
                              (bit4.be2lambda2*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda2-der.p2.lambda2))/p00 )  
  
  d2l.rho.lambda1  <- -weights*((bit1.thlambda1*p11-c.copula.theta*c.copula.lambda1)/p11+
                              (bit2.thlambda1*p10-(-c.copula.theta*(der.p1.lambda1-c.copula.lambda1)))/p10+
                              (bit3.thlambda1*p01-(-c.copula.theta*(-c.copula.lambda1)))/p01+
                              (bit4.thlambda1*p00-c.copula.theta*(c.copula.lambda1-der.p1.lambda1))/p00  )
    
  d2l.rho.lambda2  <- -weights*((bit1.thlambda2*p11-c.copula.theta*c.copula.lambda2)/p11+
                              (bit2.thlambda2*p10-(-c.copula.theta*(-c.copula.lambda2)))/p10+
                              (bit3.thlambda2*p01-(-c.copula.theta*(der.p2.lambda2-c.copula.lambda2)))/p01+
                              (bit4.thlambda2*p00-c.copula.theta*(c.copula.lambda2-der.p2.lambda2))/p00  )
  
  d2l.lambda1.lambda1  <-   -weights*((bit1.lambda1.2*p11-c.copula.lambda1^2)/p11+
                              (bit2.lambda1.2*p10-(der.p1.lambda1-c.copula.lambda1)^2)/p10+
                              (bit3.lambda1.2*p01-(-c.copula.lambda1)^2)/p01+
                              (bit4.lambda1.2*p00-(c.copula.lambda1-der.p1.lambda1)^2)/p00)
    
  d2l.lambda2.lambda2  <- -weights*((bit1.lambda2.2*p11-c.copula.lambda2^2)/p11+
                              (bit2.lambda2.2*p10-(-c.copula.lambda2)^2)/p10+
                              (bit3.lambda2.2*p01-(der.p2.lambda2-c.copula.lambda2)^2)/p01+
                              (bit4.lambda2.2*p00-(c.copula.lambda2-der.p2.lambda2)^2)/p00)
  
  d2l.lambda1.lambda2  <- -weights*((bit1.lambda1lambda2*p11-c.copula.lambda1*c.copula.lambda2)/p11+
                              (bit2.lambda1lambda2*p10-(der.p1.lambda1-c.copula.lambda1)*(-c.copula.lambda2))/p10+
                              (bit3.lambda1lambda2*p01-(-c.copula.lambda1)*(der.p2.lambda2-c.copula.lambda2))/p01+
                              (bit4.lambda1lambda2*p00-(c.copula.lambda1-der.p1.lambda1)*(c.copula.lambda2-der.p2.lambda2))/p00)


}

if(eqPL=="first"){  
                           
  d2l.be1.lambda1  <- -weights*((bit1.be1lambda1*p11-(c.copula.be1*d.n1*c.copula.lambda1))/p11+
                              (bit2.be1lambda1*p10-((1-c.copula.be1)*d.n1)*(der.p1.lambda1-c.copula.lambda1))/p10+
                              (bit3.be1lambda1*p01-(-c.copula.be1*d.n1)*(-c.copula.lambda1))/p01+
                              (bit4.be1lambda1*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda1-der.p1.lambda1))/p00 )  
  
  d2l.be2.lambda1  <-  -weights*((bit1.be2lambda1*p11-(c.copula.be2*d.n2*c.copula.lambda1))/p11+
                              (bit2.be2lambda1*p10-(-c.copula.be2*d.n2)*(der.p1.lambda1-c.copula.lambda1))/p10+
                              (bit3.be2lambda1*p01-((1-c.copula.be2)*d.n2)*(-c.copula.lambda1))/p01+
                              (bit4.be2lambda1*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda1-der.p1.lambda1))/p00   ) 
    
  d2l.rho.lambda1  <- -weights*((bit1.thlambda1*p11-c.copula.theta*c.copula.lambda1)/p11+
                              (bit2.thlambda1*p10-(-c.copula.theta*(der.p1.lambda1-c.copula.lambda1)))/p10+
                              (bit3.thlambda1*p01-(-c.copula.theta*(-c.copula.lambda1)))/p01+
                              (bit4.thlambda1*p00-c.copula.theta*(c.copula.lambda1-der.p1.lambda1))/p00  )
    
  d2l.lambda1.lambda1  <-   -weights*((bit1.lambda1.2*p11-c.copula.lambda1^2)/p11+
                              (bit2.lambda1.2*p10-(der.p1.lambda1-c.copula.lambda1)^2)/p10+
                              (bit3.lambda1.2*p01-(-c.copula.lambda1)^2)/p01+
                              (bit4.lambda1.2*p00-(c.copula.lambda1-der.p1.lambda1)^2)/p00)
    

}


if(eqPL=="second"){                             
 
  
  d2l.be1.lambda2  <-  -weights*((bit1.be1lambda2*p11-(c.copula.be1*d.n1*c.copula.lambda2))/p11+
                              (bit2.be1lambda2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.lambda2))/p10+
                              (bit3.be1lambda2*p01-(-c.copula.be1*d.n1)*(der.p2.lambda2-c.copula.lambda2))/p01+
                              (bit4.be1lambda2*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda2-der.p2.lambda2))/p00  )


    
  d2l.be2.lambda2  <- -weights*((bit1.be2lambda2*p11-(c.copula.be2*d.n2*c.copula.lambda2))/p11+
                              (bit2.be2lambda2*p10-(-c.copula.be2*d.n2)*(-c.copula.lambda2))/p10+
                              (bit3.be2lambda2*p01-((1-c.copula.be2)*d.n2)*(der.p2.lambda2-c.copula.lambda2))/p01+
                              (bit4.be2lambda2*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda2-der.p2.lambda2))/p00 )  


    
  d2l.rho.lambda2  <- -weights*((bit1.thlambda2*p11-c.copula.theta*c.copula.lambda2)/p11+
                              (bit2.thlambda2*p10-(-c.copula.theta*(-c.copula.lambda2)))/p10+
                              (bit3.thlambda2*p01-(-c.copula.theta*(der.p2.lambda2-c.copula.lambda2)))/p01+
                              (bit4.thlambda2*p00-c.copula.theta*(c.copula.lambda2-der.p2.lambda2))/p00  )


    
  d2l.lambda2.lambda2  <- -weights*((bit1.lambda2.2*p11-c.copula.lambda2^2)/p11+
                              (bit2.lambda2.2*p10-(-c.copula.lambda2)^2)/p10+
                              (bit3.lambda2.2*p01-(der.p2.lambda2-c.copula.lambda2)^2)/p01+
                              (bit4.lambda2.2*p00-(c.copula.lambda2-der.p2.lambda2)^2)/p00)
  

}


}



}
#################

if(VC$hess==FALSE && (VC$end==1 || VC$end==2)){


if(fitPL!="fixed"){



if(eqPL=="both"){  
                           
  d2l.be1.lambda1  <- -weights*(resp1*fi*(bit1.be1lambda1*p11-(c.copula.be1*d.n1*c.copula.lambda1))/p11^2+
                                resp2*se*(bit2.be1lambda1*p10-((1-c.copula.be1)*d.n1)*(der.p1.lambda1-c.copula.lambda1))/p10^2+
                                resp3*th*(bit3.be1lambda1*p01-(-c.copula.be1*d.n1)*(-c.copula.lambda1))/p01^2+
                                resp4*fo*(bit4.be1lambda1*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda1-der.p1.lambda1))/p00^2 )  
  
  d2l.be1.lambda2  <-  -weights*(resp1*fi*(bit1.be1lambda2*p11-(c.copula.be1*d.n1*c.copula.lambda2))/p11^2+
                              resp2*se*(bit2.be1lambda2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.lambda2))/p10^2+
                              resp3*th*(bit3.be1lambda2*p01-(-c.copula.be1*d.n1)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              resp4*fo*(bit4.be1lambda2*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda2-der.p2.lambda2))/p00^2  )
  
  d2l.be2.lambda1  <-  -weights*(resp1*fi*(bit1.be2lambda1*p11-(c.copula.be2*d.n2*c.copula.lambda1))/p11^2+
                              resp2*se*(bit2.be2lambda1*p10-(-c.copula.be2*d.n2)*(der.p1.lambda1-c.copula.lambda1))/p10^2+
                              resp3*th*(bit3.be2lambda1*p01-((1-c.copula.be2)*d.n2)*(-c.copula.lambda1))/p01^2+
                              resp4*fo*(bit4.be2lambda1*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda1-der.p1.lambda1))/p00^2   ) 
    
  d2l.be2.lambda2  <- -weights*(resp1*fi*(bit1.be2lambda2*p11-(c.copula.be2*d.n2*c.copula.lambda2))/p11^2+
                              resp2*se*(bit2.be2lambda2*p10-(-c.copula.be2*d.n2)*(-c.copula.lambda2))/p10^2+
                              resp3*th*(bit3.be2lambda2*p01-((1-c.copula.be2)*d.n2)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              resp4*fo*(bit4.be2lambda2*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda2-der.p2.lambda2))/p00^2 )  
  
  d2l.rho.lambda1  <- -weights*(resp1*fi*(bit1.thlambda1*p11-c.copula.theta*c.copula.lambda1)/p11^2+
                              resp2*se*(bit2.thlambda1*p10-(-c.copula.theta*(der.p1.lambda1-c.copula.lambda1)))/p10^2+
                              resp3*th*(bit3.thlambda1*p01-(-c.copula.theta*(-c.copula.lambda1)))/p01^2+
                              resp4*fo*(bit4.thlambda1*p00-c.copula.theta*(c.copula.lambda1-der.p1.lambda1))/p00^2  )
    
  d2l.rho.lambda2  <- -weights*(resp1*fi*(bit1.thlambda2*p11-c.copula.theta*c.copula.lambda2)/p11^2+
                              resp2*se*(bit2.thlambda2*p10-(-c.copula.theta*(-c.copula.lambda2)))/p10^2+
                              resp3*th*(bit3.thlambda2*p01-(-c.copula.theta*(der.p2.lambda2-c.copula.lambda2)))/p01^2+
                              resp4*fo*(bit4.thlambda2*p00-c.copula.theta*(c.copula.lambda2-der.p2.lambda2))/p00^2  )
  
  d2l.lambda1.lambda1  <-   -weights*(resp1*fi*(bit1.lambda1.2*p11-c.copula.lambda1^2)/p11^2+
                              resp2*se*(bit2.lambda1.2*p10-(der.p1.lambda1-c.copula.lambda1)^2)/p10^2+
                              resp3*th*(bit3.lambda1.2*p01-(-c.copula.lambda1)^2)/p01^2+
                              resp4*fo*(bit4.lambda1.2*p00-(c.copula.lambda1-der.p1.lambda1)^2)/p00^2)
    
  d2l.lambda2.lambda2  <- -weights*(resp1*fi*(bit1.lambda2.2*p11-c.copula.lambda2^2)/p11^2+
                              resp2*se*(bit2.lambda2.2*p10-(-c.copula.lambda2)^2)/p10^2+
                              resp3*th*(bit3.lambda2.2*p01-(der.p2.lambda2-c.copula.lambda2)^2)/p01^2+
                              resp4*fo*(bit4.lambda2.2*p00-(c.copula.lambda2-der.p2.lambda2)^2)/p00^2)
  
  d2l.lambda1.lambda2  <- -weights*(resp1*fi*(bit1.lambda1lambda2*p11-c.copula.lambda1*c.copula.lambda2)/p11^2+
                              resp2*se*(bit2.lambda1lambda2*p10-(der.p1.lambda1-c.copula.lambda1)*(-c.copula.lambda2))/p10^2+
                              resp3*th*(bit3.lambda1lambda2*p01-(-c.copula.lambda1)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              resp4*fo*(bit4.lambda1lambda2*p00-(c.copula.lambda1-der.p1.lambda1)*(c.copula.lambda2-der.p2.lambda2))/p00^2)


}

if(eqPL=="first"){  
                           
  d2l.be1.lambda1  <- -weights*(resp1*fi*(bit1.be1lambda1*p11-(c.copula.be1*d.n1*c.copula.lambda1))/p11^2+
                              resp2*se*(bit2.be1lambda1*p10-((1-c.copula.be1)*d.n1)*(der.p1.lambda1-c.copula.lambda1))/p10^2+
                              resp3*th*(bit3.be1lambda1*p01-(-c.copula.be1*d.n1)*(-c.copula.lambda1))/p01^2+
                              resp4*fo*(bit4.be1lambda1*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda1-der.p1.lambda1))/p00^2 )  
  
  d2l.be2.lambda1  <-  -weights*(resp1*fi*(bit1.be2lambda1*p11-(c.copula.be2*d.n2*c.copula.lambda1))/p11^2+
                              resp2*se*(bit2.be2lambda1*p10-(-c.copula.be2*d.n2)*(der.p1.lambda1-c.copula.lambda1))/p10^2+
                              resp3*th*(bit3.be2lambda1*p01-((1-c.copula.be2)*d.n2)*(-c.copula.lambda1))/p01^2+
                              resp4*fo*(bit4.be2lambda1*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda1-der.p1.lambda1))/p00^2   ) 
    
  d2l.rho.lambda1  <- -weights*(resp1*fi*(bit1.thlambda1*p11-c.copula.theta*c.copula.lambda1)/p11^2+
                              resp2*se*(bit2.thlambda1*p10-(-c.copula.theta*(der.p1.lambda1-c.copula.lambda1)))/p10^2+
                              resp3*th*(bit3.thlambda1*p01-(-c.copula.theta*(-c.copula.lambda1)))/p01^2+
                              resp4*fo*(bit4.thlambda1*p00-c.copula.theta*(c.copula.lambda1-der.p1.lambda1))/p00^2  )
    
  d2l.lambda1.lambda1  <-   -weights*(resp1*fi*(bit1.lambda1.2*p11-c.copula.lambda1^2)/p11^2+
                              resp2*se*(bit2.lambda1.2*p10-(der.p1.lambda1-c.copula.lambda1)^2)/p10^2+
                              resp3*th*(bit3.lambda1.2*p01-(-c.copula.lambda1)^2)/p01^2+
                              resp4*fo*(bit4.lambda1.2*p00-(c.copula.lambda1-der.p1.lambda1)^2)/p00^2)
    

}


if(eqPL=="second"){                             
 
  
  d2l.be1.lambda2  <-  -weights*(resp1*fi*(bit1.be1lambda2*p11-(c.copula.be1*d.n1*c.copula.lambda2))/p11^2+
                              resp2*se*(bit2.be1lambda2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.lambda2))/p10^2+
                              resp3*th*(bit3.be1lambda2*p01-(-c.copula.be1*d.n1)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              resp4*fo*(bit4.be1lambda2*p00-((c.copula.be1-1)*d.n1)*(c.copula.lambda2-der.p2.lambda2))/p00^2  )


    
  d2l.be2.lambda2  <- -weights*(resp1*fi*(bit1.be2lambda2*p11-(c.copula.be2*d.n2*c.copula.lambda2))/p11^2+
                              resp2*se*(bit2.be2lambda2*p10-(-c.copula.be2*d.n2)*(-c.copula.lambda2))/p10^2+
                              resp3*th*(bit3.be2lambda2*p01-((1-c.copula.be2)*d.n2)*(der.p2.lambda2-c.copula.lambda2))/p01^2+
                              resp4*fo*(bit4.be2lambda2*p00-((c.copula.be2-1)*d.n2)*(c.copula.lambda2-der.p2.lambda2))/p00^2 )  


    
  d2l.rho.lambda2  <- -weights*(resp1*fi*(bit1.thlambda2*p11-c.copula.theta*c.copula.lambda2)/p11^2+
                              resp2*se*(bit2.thlambda2*p10-(-c.copula.theta*(-c.copula.lambda2)))/p10^2+
                              resp3*th*(bit3.thlambda2*p01-(-c.copula.theta*(der.p2.lambda2-c.copula.lambda2)))/p01^2+
                              resp4*fo*(bit4.thlambda2*p00-c.copula.theta*(c.copula.lambda2-der.p2.lambda2))/p00^2  )


    
  d2l.lambda2.lambda2  <- -weights*(resp1*fi*(bit1.lambda2.2*p11-c.copula.lambda2^2)/p11^2+
                              resp2*se*(bit2.lambda2.2*p10-(-c.copula.lambda2)^2)/p10^2+
                              resp3*th*(bit3.lambda2.2*p01-(der.p2.lambda2-c.copula.lambda2)^2)/p01^2+
                              resp4*fo*(bit4.lambda2.2*p00-(c.copula.lambda2-der.p2.lambda2)^2)/p00^2)
  

}


}



}
########################







  be1.be1 <- crossprod(X1*c(d2l.be1.be1),X1)
  be2.be2 <- crossprod(X2*c(d2l.be2.be2),X2)
  be1.be2 <- crossprod(X1*c(d2l.be1.be2),X2)
  be1.rho <- t(t(rowSums(t(X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(X2*c(d2l.be2.rho)))))
  
  
  
if(fitPL!="fixed"){  
  
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

}


  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) )
            ) 

         G   <- c( colSums( c(dl.dbe1)*X1 ),
                   colSums( c(dl.dbe2)*X2 ),
                   sum( dl.drho )  ) 


if(fitPL!="fixed"){


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




}




  if(VC$extra.regI == TRUE){
  
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




    res <- -sum(l.par)
    
    if(eqPL=="both"){   add.z <- diag(c(1,1)); add.z <- add.z*c(sp.xi1,sp.xi2); nsh <- 2  } 
    if(eqPL=="first"){  add.z <- sp.xi1; nsh <- 1}
    if(eqPL=="second"){ add.z <- sp.xi2; nsh <- 1}

ps <- list()

if( ( VC$l.sp1==0 && VC$l.sp2==0 ) || VC$fp==TRUE){ lps <- length(params) - nsh 
                                                    S.h <- adiag( matrix(0,lps,lps), add.z) 
    
                                                   }else S.h <- penPL(qu.mag, sp, VC, add.z, fitPL)$S.h
         
   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
   
         if(fitPL=="fixed") G <- -G
         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
        
         H   <- H + S.h  

ps$S.h <- S.h; ps$S.h1 <- S.h1; ps$S.h2 <- S.h2  


         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps,
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
              d2l.lambda1.lambda2=d2l.lambda1.lambda2, 
              good=good, PL=PL, eqPL=eqPL, BivD=VC$BivD, p1=p1, p2=p2)      

}




     
























bprobgHsCont <- function(params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){

X1 <- X2 <- X3 <- X4 <- 1

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- etas <- NULL 

  epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)
  max.p   <- 0.9999999
  
  
if(is.null(VC$X3)){  
  sigma2.st <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
  teta.st   <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
} 

if(!is.null(VC$X3)){  
  sigma2.st <- etas <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  teta.st   <- etad <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
}  
  
  sigma2 <- exp(sigma2.st) + epsilon

  cy <- respvec$cy  
  y2 <- respvec$y2 
  y1 <- respvec$y1   
  weights <- VC$weights
  
 dHs  <- distrHs(y2, eta2, sigma2, sigma2.st, margin2=VC$margins[2], naive = FALSE)
  
 pdf2                         <- dHs$pdf2
 p2                           <- dHs$p2 
 derpdf2.dereta2              <- dHs$derpdf2.dereta2 
 derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st 
 derp2.dersigma.st            <- dHs$derp2.dersigma.st
 derp2.dereta2                <- dHs$derp2.dereta2
 der2p2.dereta2eta2           <- dHs$der2p2.dereta2eta2 
 der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
 der2p2.dersigma2.st2         <- dHs$der2p2.dersigma2.st2
 der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2
 der2p2.dereta2dersigma2.st   <- dHs$der2p2.dereta2dersigma2.st            
 der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  
  
  
  p1 <- pnorm(-eta1)
  p1 <- ifelse(p1 > max.p, max.p, p1) 
  p1 <- ifelse(p1 < epsilon, epsilon, p1) 

########################################################################################################  
  
  if( VC$BivD %in% c("N") ) teta.st <- ifelse(abs(teta.st) > 8.75, sign(teta.st)*8.75, teta.st )  
  
  if( VC$BivD %in%  c("C0","C180","C90","C270","J0","J180","J90","J270","G0","G180","G90","G270") ) {
 
  teta.st <- ifelse( teta.st > 20, 20, teta.st )  # 709
  teta.st <- ifelse( teta.st < -17, -17, teta.st )   # -20
  
  }
  
 
    if(VC$BivD %in% c("N")      )     teta <- tanh(teta.st)                        
    if(VC$BivD=="F")                  teta <- teta.st + epsilon
    if(VC$BivD %in% c("C0", "C180") ) teta <-    exp(teta.st)                      
    if(VC$BivD %in% c("C90","C270") ) teta <- -( exp(teta.st) )                     
    if(VC$BivD %in% c("J0", "J180","G0", "G180") ) teta <-    exp(teta.st) + 1     
    if(VC$BivD %in% c("J90","J270","G90","G270") ) teta <- -( exp(teta.st) + 1 )   
    

########################################################################################################

dH <- copgHs(p1,p2,eta1=NULL,eta2=NULL,teta,teta.st,VC$BivD)

h <- dH$c.copula.be2  
h <- ifelse(h > max.p, max.p, h) 
h <- ifelse(h < epsilon  , epsilon  , h) 


l.par <- weights*( cy*log(h) + y1*log(1 - h) + log(pdf2) )
  
########################################################################################################
 
  c.copula2.be2    <- dH$c.copula2.be2 
  c.copula2.be1be2 <- dH$c.copula2.be1be2   
  c.copula2.be2th  <- dH$c.copula2.be2th  
  derp1.dereta1    <- -dnorm(-eta1) 
  derh.dereta1     <- c.copula2.be1be2 * derp1.dereta1
  derh.dereta2     <- c.copula2.be2 * derp2.dereta2

  dl.dbe1      <- weights*( derh.dereta1 *(cy/h -y1/(1-h)) ) 
  dl.dbe2      <- weights*( derh.dereta2 * (cy/h -y1/(1-h)) + derpdf2.dereta2/pdf2 )
  dl.dsigma.st <- weights*( c.copula2.be2 * derp2.dersigma.st *(cy/h -y1/(1-h)) + derpdf2.dersigma2.st/pdf2 )
  dl.dteta.st  <- weights*( c.copula2.be2th*(cy/h -y1/(1-h)) )                     
 
 
######################################################################################################## 
 
  BITS <- copgHsCont(p1, p2, teta, teta.st, VC)
  
  der2h.derp2p2              <- BITS$der2h.derp2p2 
  der2h.derteta.teta.st      <- BITS$der2h.derteta.teta.st  
  derteta.derteta.st         <- BITS$derteta.derteta.st 
  der2teta.derteta.stteta.st <- BITS$der2teta.derteta.stteta.st  
  der2h.derp1p2              <- BITS$der2h.derp1p2  
  der2h.derp1teta            <- BITS$der2h.derp1teta                                     
  der2h.derp2teta            <- BITS$der2h.derp2teta  
  der2h.derp1p1              <- BITS$der2h.derp1p1
  
 
  der2p1.dereta1eta1 <- eta1 * dnorm(-eta1)      
                    
  der2h.dereta2.dereta2         <- der2h.derp2p2*derp2.dereta2^2 + c.copula2.be2*der2p2.dereta2eta2                                        
  der2h.derteta.st2             <- der2h.derteta.teta.st*derteta.derteta.st^2 + c.copula2.be2th/derteta.derteta.st* der2teta.derteta.stteta.st                                                                                        
  der2h.derp2dersigma2.st       <- der2h.derp2p2*derp2.dersigma.st 
  der2h.dersigma2.st2           <- der2h.derp2dersigma2.st*derp2.dersigma.st + c.copula2.be2*der2p2.dersigma2.st2
  derh.dersigma2.st             <- c.copula2.be2 * derp2.dersigma.st
  der2h.dereta1.dereta2         <- der2h.derp1p2*derp1.dereta1*derp2.dereta2                                                                   
  der2h.dereta1.derteta.st      <- der2h.derp1teta*derp1.dereta1*derteta.derteta.st  
  der2h.dereta1.dersigma2.st    <- der2h.derp1p2 * derp2.dersigma.st*derp1.dereta1   
  der2h.dereta2.derteta.st      <- der2h.derp2teta*derp2.dereta2*derteta.derteta.st 
  der2h.derteta.st.dersigma2.st <- der2h.derp2teta* derteta.derteta.st*derp2.dersigma.st  
  der2h.dereta2.dersigma2.st    <- der2h.derp2dersigma2.st*derp2.dereta2 + c.copula2.be2*der2p2.dereta2dersigma2.st                                        
  der2h.dereta1.dereta1         <- der2h.derp1p1*derp1.dereta1^2 + c.copula2.be1be2*der2p1.dereta1eta1      
  
  
  d2l.be1.be1      <- -weights*(der2h.dereta1.dereta1 *(cy/h -y1/(1-h)) - derh.dereta1^2 * (cy/h^2 +y1/(1-h)^2) )
  d2l.be2.be2      <- -weights*(der2h.dereta2.dereta2 *(cy/h -y1/(1-h)) - derh.dereta2^2 * (cy/h^2 +y1/(1-h)^2) + (der2pdf2.dereta2*pdf2-(derpdf2.dereta2)^2)/(pdf2)^2 )
  d2l.rho.rho      <- -weights*(der2h.derteta.st2 *(cy/h -y1/(1-h)) - c.copula2.be2th^2 * (cy/h^2 +y1/(1-h)^2) )
  d2l.sigma.sigma  <- -weights*(der2h.dersigma2.st2 *(cy/h -y1/(1-h)) - derh.dersigma2.st^2 * (cy/h^2 +y1/(1-h)^2) + (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/(pdf2)^2 )

  d2l.be1.be2      <- -weights*(der2h.dereta1.dereta2 *(cy/h -y1/(1-h)) - derh.dereta1*derh.dereta2 * (cy/h^2 +y1/(1-h)^2)  )
  d2l.be1.rho      <- -weights*(der2h.dereta1.derteta.st *(cy/h -y1/(1-h)) - derh.dereta1 * c.copula2.be2th* (cy/h^2 +y1/(1-h)^2) )
  d2l.be1.sigma    <- -weights*(der2h.dereta1.dersigma2.st *(cy/h -y1/(1-h)) - derh.dereta1 * derh.dersigma2.st* (cy/h^2 +y1/(1-h)^2) )
  
  d2l.be2.rho      <- -weights*(der2h.dereta2.derteta.st *(cy/h -y1/(1-h)) - derh.dereta2*c.copula2.be2th * (cy/h^2 +y1/(1-h)^2)  )
  d2l.be2.sigma    <- -weights*(der2h.dereta2.dersigma2.st *(cy/h -y1/(1-h)) - derh.dereta2*derh.dersigma2.st * (cy/h^2 +y1/(1-h)^2) + (der2pdf2.dereta2dersigma2.st*pdf2-(derpdf2.dereta2*derpdf2.dersigma2.st))/(pdf2)^2 ) 
  
  d2l.rho.sigma    <- -weights*(der2h.derteta.st.dersigma2.st *(cy/h -y1/(1-h)) - derh.dersigma2.st * c.copula2.be2th* (cy/h^2 +y1/(1-h)^2) )
  

if( is.null(VC$X3) ){


  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho   <- t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be1.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma))))) 
  be2.rho   <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma))))) 

  H <- rbind( cbind( be1.be1    ,   be1.be2    ,   be1.sigma,            be1.rho  ), 
              cbind( t(be1.be2) ,   be2.be2    ,   be2.sigma,            be2.rho  ), 
              cbind( t(be1.sigma),  t(be2.sigma),  sum(d2l.sigma.sigma), sum(d2l.rho.sigma) ),
              cbind( t(be1.rho) ,   t(be2.rho),    sum(d2l.rho.sigma),   sum(d2l.rho.rho)   ) 
              ) 
         
  G   <- -c( colSums( c(dl.dbe1)*VC$X1 ) ,
             colSums( c(dl.dbe2)*VC$X2 ) ,
             sum( dl.dsigma.st ),
             sum( dl.dteta.st ) )
    
}




if( !is.null(VC$X3) ){

  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  
  be1.rho   <- crossprod(VC$X1*c(d2l.be1.rho),  VC$X4)                                     
  be1.sigma <- crossprod(VC$X1*c(d2l.be1.sigma),VC$X3)                                   
  be2.rho   <- crossprod(VC$X2*c(d2l.be2.rho),  VC$X4)                                     
  be2.sigma <- crossprod(VC$X2*c(d2l.be2.sigma),VC$X3)  
  
  sigma.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)   
  sigma.rho   <- crossprod(VC$X3*c(d2l.rho.sigma),VC$X4)  
  rho.rho     <- crossprod(VC$X4*c(d2l.rho.rho),    VC$X4)    
  
  

  H <- rbind( cbind( be1.be1     ,  be1.be2     ,  be1.sigma   , be1.rho   ), 
              cbind( t(be1.be2)  ,  be2.be2     ,  be2.sigma   , be2.rho   ), 
              cbind( t(be1.sigma),  t(be2.sigma),  sigma.sigma , sigma.rho ),
              cbind( t(be1.rho)  ,  t(be2.rho)  ,  t(sigma.rho), rho.rho   ) 
             )  
            
   
  G   <- -c( colSums(      c(dl.dbe1)*VC$X1 ) ,
             colSums(      c(dl.dbe2)*VC$X2 ) ,
             colSums( c(dl.dsigma.st)*VC$X3 ) ,
             colSums(  c(dl.dteta.st)*VC$X4 ) )   
   
    
}




      res <- -sum(l.par)




if( ( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0 && VC$l.sp4==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)

 
if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  
  
rm(X1,X2,X3,X4)  
  
  
  

         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps, etas = etas,
              eta1=eta1, eta2=eta2, etad=etad,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.dsigma.st = dl.dsigma.st, dl.dteta.st = dl.dteta.st,
d2l.be1.be1      =d2l.be1.be1 ,   
d2l.be2.be2    	 =d2l.be2.be2 ,   
d2l.rho.rho    	 =d2l.rho.rho ,   
d2l.sigma.sigma	 =d2l.sigma.sigma,
d2l.be1.be2    	 =d2l.be1.be2    ,
d2l.be1.rho    	 =d2l.be1.rho    ,
d2l.be1.sigma  	 =d2l.be1.sigma  ,
d2l.be2.rho    	 =d2l.be2.rho    ,
d2l.be2.sigma  	 =d2l.be2.sigma  ,
d2l.rho.sigma  	 =d2l.rho.sigma  ,
              BivD=VC$BivD, p1=1-p1, p2=p2, theta.star = teta.st)      

}




     
























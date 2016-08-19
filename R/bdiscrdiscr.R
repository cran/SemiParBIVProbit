bdiscrdiscr <- function(params, respvec, VC, ps, AT = FALSE){

    eta1 <- VC$X1%*%params[1:VC$X1.d2]
    eta2 <- VC$X2%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
    etad <- etas1 <- etas2 <- l.ln <- NULL 
  
    epsilon <- 0.0000001 
    max.p   <- 0.9999999
    
    
  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
  } 
  
  
  if(!is.null(VC$X3)){  
    sigma21.st <- etas1 <- VC$X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- VC$X4%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    teta.st    <- etad  <- VC$X5%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
  }  
  
  
##################
## Transformations
##################
  
    sstr1 <- esp.tr(sigma21.st, VC$margins[1])  
    sstr2 <- esp.tr(sigma22.st, VC$margins[2])  
  
    sigma21.st <- sstr1$vrb.st 
    sigma22.st <- sstr2$vrb.st 
    
    sigma21    <- sstr1$vrb 
    sigma22    <- sstr2$vrb 
    
    eta1 <- eta.tr(eta1, VC$margins[1])
    eta2 <- eta.tr(eta2, VC$margins[2])
    
resT    <- teta.tr(VC, teta.st)
teta.st <- resT$teta.st
teta    <- resT$teta
    
##################
##################

  dHs1 <- distrHsDiscr(respvec$y1, eta1, sigma21, sigma21.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = FALSE, y2m = VC$y1m)
  dHs2 <- distrHsDiscr(respvec$y2, eta2, sigma22, sigma22.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = FALSE, y2m = VC$y2m)

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
  
  C11 <- BiCDF(p1, p2, VC$nC, teta)
  C01 <- BiCDF(mm(p1-pdf1), p2, VC$nC, teta)
  C10 <- BiCDF(p1, mm(p2-pdf2), VC$nC, teta)
  C00 <- BiCDF(mm(p1-pdf1), mm(p2-pdf2), VC$nC, teta)

  E <- C11 - C01 - C10 + C00 
  E <- ifelse(E < epsilon, epsilon, E)  


  l.par <- VC$weights*log(E)
  
  
##################

  dHC11 <- copgHs(p1, p2,                   eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD)
  dHC01 <- copgHs(mm(p1-pdf1), p2,          eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD)
  dHC10 <- copgHs(p1, mm(p2-pdf2),          eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD)
  dHC00 <- copgHs(mm(p1-pdf1), mm(p2-pdf2), eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD)
 
  derC11.derp1 <- dHC11$c.copula.be1  
  derC01.derp1 <- dHC01$c.copula.be1  
  derC10.derp1 <- dHC10$c.copula.be1  
  derC00.derp1 <- dHC00$c.copula.be1 
  
  derC11.derp2 <- dHC11$c.copula.be2  
  derC01.derp2 <- dHC01$c.copula.be2  
  derC10.derp2 <- dHC10$c.copula.be2  
  derC00.derp2 <- dHC00$c.copula.be2   
  
  derC11.derthet <- dHC11$c.copula.thet 
  derC01.derthet <- dHC01$c.copula.thet  
  derC10.derthet <- dHC10$c.copula.thet  
  derC00.derthet <- dHC00$c.copula.thet  
  
  derteta.derteta.st <- dHC11$derteta.derteta.st
  
  derpdf1.dereta1    <- dHs1$derpdf2.dereta2 
  derp1.dereta1      <- dHs1$derp2.dereta2
  derp1m1.dereta1    <- derp1.dereta1 - derpdf1.dereta1
  
  derpdf1.dersigma21.st  <- dHs1$derpdf2.dersigma2.st  
  derp1.dersigma21.st    <- dHs1$derp2.dersigma.st 
  derp1m1.dersigma21.st  <- derp1.dersigma21.st - derpdf1.dersigma21.st    
  
  derpdf2.dereta2    <- dHs2$derpdf2.dereta2 
  derp2.dereta2      <- dHs2$derp2.dereta2
  derp2m1.dereta2    <- derp2.dereta2 - derpdf2.dereta2
  
  derpdf2.dersigma22.st  <- dHs2$derpdf2.dersigma2.st  
  derp2.dersigma22.st    <- dHs2$derp2.dersigma.st 
  derp2m1.dersigma22.st  <- derp2.dersigma22.st - derpdf2.dersigma22.st      

  fE1  <- (derC11.derp1 - derC10.derp1)*derp1.dereta1 - (derC01.derp1 - derC00.derp1)*derp1m1.dereta1    
  fE2  <- (derC11.derp2 - derC01.derp2)*derp2.dereta2 - (derC10.derp2 - derC00.derp2)*derp2m1.dereta2   
  fEt  <- (derC11.derthet - derC01.derthet - derC10.derthet + derC00.derthet)*derteta.derteta.st   
  fE1s <- (derC11.derp1 - derC10.derp1)*derp1.dersigma21.st - (derC01.derp1 - derC00.derp1)*derp1m1.dersigma21.st    
  fE2s <- (derC11.derp2 - derC01.derp2)*derp2.dersigma22.st - (derC10.derp2 - derC00.derp2)*derp2m1.dersigma22.st   

  dl.dbe1        <- VC$weights*fE1/E   
  dl.dbe2        <- VC$weights*fE2/E
  dl.dsigma21.st <- VC$weights*fE1s/E
  dl.dsigma22.st <- VC$weights*fE2s/E
  dl.dteta.st    <- VC$weights*fEt/E

##################

  der2C11.derp1p1 <- dHC11$c.copula2.be1
  der2C01.derp1p1 <- dHC01$c.copula2.be1  
  der2C10.derp1p1 <- dHC10$c.copula2.be1  
  der2C00.derp1p1 <- dHC00$c.copula2.be1 
  
  der2C11.derp2p2 <- dHC11$c.copula2.be2
  der2C01.derp2p2 <- dHC01$c.copula2.be2  
  der2C10.derp2p2 <- dHC10$c.copula2.be2  
  der2C00.derp2p2 <- dHC00$c.copula2.be2  
  
  der2C11.derp1p2 <- dHC11$c.copula2.be1be2
  der2C01.derp1p2 <- dHC01$c.copula2.be1be2 
  der2C10.derp1p2 <- dHC10$c.copula2.be1be2  
  der2C00.derp1p2 <- dHC00$c.copula2.be1be2
 
  der2C11.derp1t  <- dHC11$c.copula2.be1t
  der2C01.derp1t  <- dHC01$c.copula2.be1t 
  der2C10.derp1t  <- dHC10$c.copula2.be1t  
  der2C00.derp1t  <- dHC00$c.copula2.be1t 
  
  der2C11.derp2t  <- dHC11$c.copula2.be2t
  der2C01.derp2t  <- dHC01$c.copula2.be2t 
  der2C10.derp2t  <- dHC10$c.copula2.be2t  
  der2C00.derp2t  <- dHC00$c.copula2.be2t   
  
  der2C11.derthet2 <- dHC11$bit1.th2ATE
  der2C01.derthet2 <- dHC01$bit1.th2ATE  
  der2C10.derthet2 <- dHC10$bit1.th2ATE  
  der2C00.derthet2 <- dHC00$bit1.th2ATE   

  der2teta.derteta.stteta.st <- dHC11$der2teta.derteta.stteta.st 

  der2pdf1.dereta1     <- dHs1$der2pdf2.dereta2 
  der2p1.dereta1eta1   <- dHs1$der2p2.dereta2eta2
  der2p1m1.dereta1eta1 <- der2p1.dereta1eta1 - der2pdf1.dereta1

  der2pdf2.dereta2     <- dHs2$der2pdf2.dereta2
  der2p2.dereta2eta2   <- dHs2$der2p2.dereta2eta2
  der2p2m1.dereta2eta2 <- der2p2.dereta2eta2 - der2pdf2.dereta2

  der2pdf1.dersigma21.st2 <- dHs1$der2pdf2.dersigma2.st2
  der2p1.dersigma21.st2   <- dHs1$der2p2.dersigma2.st2
  der2p1m1.dersigma21.st2 <- der2p1.dersigma21.st2 - der2pdf1.dersigma21.st2  
  
  der2pdf2.dersigma22.st2 <- dHs2$der2pdf2.dersigma2.st2
  der2p2.dersigma22.st2   <- dHs2$der2p2.dersigma2.st2
  der2p2m1.dersigma22.st2 <- der2p2.dersigma22.st2 - der2pdf2.dersigma22.st2  

  der2pdf1.dereta1dersigma21.st <- dHs1$der2pdf2.dereta2dersigma2.st 
  der2p1.dereta1dersigma21.st   <- dHs1$der2p2.dereta2dersigma2.st
  der2p1m1.dereta1dersigma21.st <- der2p1.dereta1dersigma21.st - der2pdf1.dereta1dersigma21.st
  
  der2pdf2.dereta2dersigma22.st <- dHs2$der2pdf2.dereta2dersigma2.st
  der2p2.dereta2dersigma22.st   <- dHs2$der2p2.dereta2dersigma2.st
  der2p2m1.dereta2dersigma22.st <- der2p2.dereta2dersigma22.st - der2pdf2.dereta2dersigma22.st

d2l.be1.be1 <- -VC$weights*( (E*(der2C11.derp1p1*derp1.dereta1^2 + derC11.derp1*der2p1.dereta1eta1 - 
(der2C01.derp1p1*derp1m1.dereta1^2 + derC01.derp1*der2p1m1.dereta1eta1) - 
(der2C10.derp1p1*derp1.dereta1^2 + derC10.derp1*der2p1.dereta1eta1) + 
(der2C00.derp1p1*derp1m1.dereta1^2 + derC00.derp1*der2p1m1.dereta1eta1)) - fE1^2)/E^2  )

d2l.be2.be2 <- -VC$weights*((E*(der2C11.derp2p2*derp2.dereta2^2 + derC11.derp2*der2p2.dereta2eta2 - 
(der2C01.derp2p2*derp2.dereta2^2 + derC01.derp2*der2p2.dereta2eta2) - 
(der2C10.derp2p2*derp2m1.dereta2^2 + derC10.derp2*der2p2m1.dereta2eta2) + 
(der2C00.derp2p2*derp2m1.dereta2^2 + derC00.derp2*der2p2m1.dereta2eta2)) - fE2^2)/E^2  )

d2l.sigma21.sigma21 <- -VC$weights*( (E*(der2C11.derp1p1*derp1.dersigma21.st^2 + derC11.derp1*der2p1.dersigma21.st2 - 
(der2C01.derp1p1*derp1m1.dersigma21.st^2 + derC01.derp1*der2p1m1.dersigma21.st2) - 
(der2C10.derp1p1*derp1.dersigma21.st^2 + derC10.derp1*der2p1.dersigma21.st2) + 
(der2C00.derp1p1*derp1m1.dersigma21.st^2 + derC00.derp1*der2p1m1.dersigma21.st2)) - fE1s^2)/E^2  )


d2l.sigma22.sigma22 <- -VC$weights*((E*(der2C11.derp2p2*derp2.dersigma22.st^2 + derC11.derp2*der2p2.dersigma22.st2 - 
(der2C01.derp2p2*derp2.dersigma22.st^2 + derC01.derp2*der2p2.dersigma22.st2) - 
(der2C10.derp2p2*derp2m1.dersigma22.st^2 + derC10.derp2*der2p2m1.dersigma22.st2) + 
(der2C00.derp2p2*derp2m1.dersigma22.st^2 + derC00.derp2*der2p2m1.dersigma22.st2)) - fE2s^2)/E^2  )

d2l.rho.rho <- -VC$weights*((E*(der2C11.derthet2*derteta.derteta.st^2 + derC11.derthet*der2teta.derteta.stteta.st - 
(der2C01.derthet2*derteta.derteta.st^2 + derC01.derthet*der2teta.derteta.stteta.st) - 
(der2C10.derthet2*derteta.derteta.st^2 + derC10.derthet*der2teta.derteta.stteta.st) + 
(der2C00.derthet2*derteta.derteta.st^2 + derC00.derthet*der2teta.derteta.stteta.st)) - fEt^2)/E^2  )                                   

d2l.be1.be2 <- -VC$weights*( (E*(der2C11.derp1p2*derp1.dereta1*derp2.dereta2  - 
der2C01.derp1p2*derp1m1.dereta1*derp2.dereta2 - 
der2C10.derp1p2*derp1.dereta1*derp2m1.dereta2 + 
der2C00.derp1p2*derp1m1.dereta1*derp2m1.dereta2) - fE1*fE2)/E^2  )

d2l.be1.sigma22 <- -VC$weights*((E*(der2C11.derp1p2*derp1.dereta1*derp2.dersigma22.st  - 
der2C01.derp1p2*derp1m1.dereta1*derp2.dersigma22.st - 
der2C10.derp1p2*derp1.dereta1*derp2m1.dersigma22.st + 
der2C00.derp1p2*derp1m1.dereta1*derp2m1.dersigma22.st) - fE1*fE2s)/E^2  )

d2l.be1.rho <- -VC$weights*((E*(der2C11.derp1t*derp1.dereta1*derteta.derteta.st  - 
der2C01.derp1t*derp1m1.dereta1*derteta.derteta.st - 
der2C10.derp1t*derp1.dereta1*derteta.derteta.st + 
der2C00.derp1t*derp1m1.dereta1*derteta.derteta.st) - fE1*fEt)/E^2  )

d2l.be1.sigma21 <- -VC$weights*((E*(der2C11.derp1p1*derp1.dereta1*derp1.dersigma21.st + derC11.derp1*der2p1.dereta1dersigma21.st - 
(der2C01.derp1p1*derp1m1.dereta1*derp1m1.dersigma21.st + derC01.derp1*der2p1m1.dereta1dersigma21.st) - 
(der2C10.derp1p1*derp1.dereta1*derp1.dersigma21.st + derC10.derp1*der2p1.dereta1dersigma21.st) + 
(der2C00.derp1p1*derp1m1.dereta1*derp1m1.dersigma21.st + derC00.derp1*der2p1m1.dereta1dersigma21.st)) - fE1*fE1s)/E^2  )


d2l.be2.sigma22 <- -VC$weights*((E*(der2C11.derp2p2*derp2.dereta2*derp2.dersigma22.st + derC11.derp2*der2p2.dereta2dersigma22.st - 
(der2C01.derp2p2*derp2.dereta2*derp2.dersigma22.st + derC01.derp2*der2p2.dereta2dersigma22.st) - 
(der2C10.derp2p2*derp2m1.dereta2*derp2m1.dersigma22.st + derC10.derp2*der2p2m1.dereta2dersigma22.st) + 
(der2C00.derp2p2*derp2m1.dereta2*derp2m1.dersigma22.st + derC00.derp2*der2p2m1.dereta2dersigma22.st)) - fE2*fE2s)/E^2
  )


d2l.be2.sigma21 <- -VC$weights*((E*(der2C11.derp1p2*derp1.dersigma21.st*derp2.dereta2  - 
der2C01.derp1p2*derp1m1.dersigma21.st*derp2.dereta2 - 
der2C10.derp1p2*derp1.dersigma21.st*derp2m1.dereta2 + 
der2C00.derp1p2*derp1m1.dersigma21.st*derp2m1.dereta2) - fE1s*fE2)/E^2  )

d2l.be2.rho <- -VC$weights*((E*(der2C11.derp2t*derp2.dereta2*derteta.derteta.st  - 
der2C01.derp2t*derp2.dereta2*derteta.derteta.st - 
der2C10.derp2t*derp2m1.dereta2*derteta.derteta.st + 
der2C00.derp2t*derp2m1.dereta2*derteta.derteta.st) - fE2*fEt)/E^2  )

d2l.rho.sigma21 <- -VC$weights*((E*(der2C11.derp1t*derp1.dersigma21.st*derteta.derteta.st  - 
der2C01.derp1t*derp1m1.dersigma21.st*derteta.derteta.st - 
der2C10.derp1t*derp1.dersigma21.st*derteta.derteta.st + 
der2C00.derp1t*derp1m1.dersigma21.st*derteta.derteta.st) - fE1s*fEt)/E^2  )

d2l.rho.sigma22 <- -VC$weights*((E*(der2C11.derp2t*derp2.dersigma22.st*derteta.derteta.st  - 
der2C01.derp2t*derp2.dersigma22.st*derteta.derteta.st - 
der2C10.derp2t*derp2m1.dersigma22.st*derteta.derteta.st + 
der2C00.derp2t*derp2m1.dersigma22.st*derteta.derteta.st) - fE2s*fEt)/E^2  )


d2l.sigma21.sigma22 <- -VC$weights*((E*(der2C11.derp1p2*derp1.dersigma21.st*derp2.dersigma22.st  - 
der2C01.derp1p2*derp1m1.dersigma21.st*derp2.dersigma22.st - 
der2C10.derp1p2*derp1.dersigma21.st*derp2m1.dersigma22.st + 
der2C00.derp1p2*derp1m1.dersigma21.st*derp2m1.dersigma22.st) - fE1s*fE2s)/E^2  )




    


if( is.null(VC$X3) ){



  G   <- -c(colSums( c(dl.dbe1)*VC$X1 ),
            colSums( c(dl.dbe2)*VC$X2 ),
            sum( dl.dsigma21.st ),
            sum( dl.dsigma22.st ),
            sum( dl.dteta.st )            )



  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <-   t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be1.sigma21 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma21))))) 
  be1.sigma22 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma22))))) 
  be2.rho <-   t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma21 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma21))))) 
  be2.sigma22 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma22)))))

  H <- rbind( cbind( be1.be1    ,   be1.be2    ,   be1.sigma21,    be1.sigma22,       be1.rho  ), 
              cbind( t(be1.be2) ,   be2.be2    ,   be2.sigma21,    be2.sigma22,        be2.rho  ), 
              cbind( t(be1.sigma21) , t(be2.sigma21) , sum(d2l.sigma21.sigma21), sum(d2l.sigma21.sigma22), sum(d2l.rho.sigma21)),
              cbind( t(be1.sigma22) , t(be2.sigma22) , sum(d2l.sigma21.sigma22), sum(d2l.sigma22.sigma22), sum(d2l.rho.sigma22)),
              cbind( t(be1.rho) ,   t(be2.rho) ,   sum(d2l.rho.sigma21), sum(d2l.rho.sigma22), sum(d2l.rho.rho) ) 
              
              ) 

}



if( !is.null(VC$X3) ){



G   <- -c( colSums(       c(dl.dbe1)*VC$X1 ) ,
           colSums(       c(dl.dbe2)*VC$X2 ) ,
           colSums(c(dl.dsigma21.st)*VC$X3 ) ,
           colSums(c(dl.dsigma22.st)*VC$X4 ) ,
           colSums(   c(dl.dteta.st)*VC$X5 )  )  

                 
    be1.be1         <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
    be2.be2         <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
    be1.be2         <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
    be1.rho         <- crossprod(VC$X1*c(d2l.be1.rho),VC$X5)    
    be2.rho         <- crossprod(VC$X2*c(d2l.be2.rho),VC$X5)   
    be1.sigma21     <- crossprod(VC$X1*c(d2l.be1.sigma21),VC$X3)      
    be1.sigma22     <- crossprod(VC$X1*c(d2l.be1.sigma22),VC$X4)      
    be2.sigma21     <- crossprod(VC$X2*c(d2l.be2.sigma21),VC$X3)   
    be2.sigma22     <- crossprod(VC$X2*c(d2l.be2.sigma22),VC$X4)   
    sigma21.sigma21 <- crossprod(VC$X3*c(d2l.sigma21.sigma21),VC$X3)     
    sigma21.sigma22 <- crossprod(VC$X3*c(d2l.sigma21.sigma22),VC$X4)
    rho.sigma21     <- crossprod(VC$X3*c(d2l.rho.sigma21),VC$X5)
    sigma22.sigma22 <- crossprod(VC$X4*c(d2l.sigma22.sigma22),VC$X4)
    rho.sigma22     <- crossprod(VC$X4*c(d2l.rho.sigma22),VC$X5)   
    rho.rho         <- crossprod(VC$X5*c(d2l.rho.rho),VC$X5)    
    
    
    H <- rbind( cbind( be1.be1        ,   be1.be2      ,   be1.sigma21,      be1.sigma22,       be1.rho    ), 
                cbind( t(be1.be2)     ,   be2.be2      ,   be2.sigma21,      be2.sigma22,       be2.rho    ), 
                cbind( t(be1.sigma21) , t(be2.sigma21) ,   sigma21.sigma21,  sigma21.sigma22,   rho.sigma21),
                cbind( t(be1.sigma22) , t(be2.sigma22) , t(sigma21.sigma22), sigma22.sigma22,   rho.sigma22),
                cbind( t(be1.rho)     ,   t(be2.rho)   ,   t(rho.sigma21),   t(rho.sigma22),    rho.rho    ) ) 
                
                 

}



res <- -sum(l.par)



 
if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
  S.h  <- ps$S.h  
  
  if( length(S.h) != 1){
  
  S.h1 <- 0.5*crossprod(params,S.h)%*%params
  S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h   
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  
  
  

  
  
  
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, 
              eta1=eta1, eta2=eta2, etad=etad, etas1 = etas1, etas2 = etas2, 
              BivD=VC$BivD, p1 = p1, p2 = p2,
              dl.dbe1          =dl.dbe1,       
              dl.dbe2          =dl.dbe2,       
              dl.dsigma21.st   =dl.dsigma21.st,
              dl.dsigma22.st   =dl.dsigma22.st,
              dl.dteta.st      =dl.dteta.st) 
              
              

  }
  

  

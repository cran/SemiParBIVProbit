bprobgHsBB <- function(params, BivD, nC, nu, sp.xi1, sp.xi2, PL, eqPL, H.n, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, weights=weights, X1.d2, X2.d2, pPen1=NULL, pPen2=NULL, sp=NULL, qu.mag=NULL, gp1, gp2, fp, l.sp1, l.sp2){

  eta1 <- X1%*%params[1:X1.d2]
  eta2 <- X2%*%params[(X1.d2+1):(X1.d2+X2.d2)]

  p1 <- pnorm(eta1); d.n1 <- dnorm(eta1)
  p2 <- pnorm(eta2); d.n2 <- dnorm(eta2) 

  teta.st  <- params[(X1.d2+X2.d2+1)]
  delta.st <- params[(X1.d2+X2.d2+2)]
  epsilon <- .Machine$double.eps*10^6

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

    if(BivD %in% c("BB1.0", "BB1.180")){teta <-   exp(teta.st) + epsilon;  delta <-   exp(delta.st) + 1}
    if(BivD %in% c("BB1.90","BB1.270")){teta <- -(exp(teta.st) + epsilon); delta <- -(exp(delta.st) + 1)}

    if(BivD %in% c("BB6.0", "BB6.180")){teta <-   exp(teta.st) + 1; delta <-   exp(delta.st) + 1}
    if(BivD %in% c("BB6.90","BB6.270")){teta <- -(exp(teta.st) + 1);delta <- -(exp(delta.st) + 1)}

    if(BivD %in% c("BB7.0", "BB7.180")){teta <-   exp(teta.st) + 1; delta <-   exp(delta.st) + epsilon}
    if(BivD %in% c("BB7.90","BB7.270")){teta <- -(exp(teta.st) + 1);delta <- -(exp(delta.st) + epsilon)}

    if(BivD %in% c("BB8.0", "BB8.180")){teta <-   exp(teta.st) + 1; delta <-  pnorm(delta.st); if(delta==0) delta <- epsilon}
    if(BivD %in% c("BB8.90","BB8.270")){teta <- -(exp(teta.st) + 1);delta <- -pnorm(delta.st); if(delta==0) delta <- -epsilon}
    
  C.copula <- BiCopCDF(p1,p2, nC, par=teta, par2=delta)
########################################################################################################


  p11 <- pmax( C.copula, 1000*.Machine$double.eps )
  p10 <- pmax( p1 - p11, 1000*.Machine$double.eps )
  p01 <- pmax( p2 - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )


  l.par <- weights*( y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) )

dH <- copgHsBB(p1,p2,teta,teta.st,delta,delta.st,BivD)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
c.copula.theta <- dH$c.copula.theta 
c.copula.delta <- dH$c.copula.delta 
  
c.copula2.be1 <- dH$c.copula2.be1   
c.copula2.be2 <- dH$c.copula2.be2 
c.copula2.be1be2 <- dH$c.copula2.be1be2
c.copula2.be1th <- dH$c.copula2.be1th 
c.copula2.be2th <- dH$c.copula2.be2th
c.copula2.be1del <- dH$c.copula2.be1del
c.copula2.be2del <- dH$c.copula2.be2del
bit1.thdel <- dH$bit1.thdel
bit1.th2 <- dH$bit1.th2
bit1.del2 <- dH$bit1.del2


bit1.b1b1 <- c.copula2.be1*(d.n1)^2-c.copula.be1*d.n1*eta1
bit2.b1b1 <- -d.n1*eta1-bit1.b1b1
bit3.b1b1 <- -bit1.b1b1
bit4.b1b1 <- -bit2.b1b1

bit1.b2b2 <- c.copula2.be2*(d.n2)^2-c.copula.be2*d.n2*eta2
bit2.b2b2 <- -bit1.b2b2
bit3.b2b2 <- -d.n2*eta2-bit1.b2b2
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

bit2.del2 <- -bit1.del2 
bit3.del2 <- -bit1.del2 
bit4.del2 <- bit1.del2 

bit1.b1del <- c.copula2.be1del*d.n1
bit2.b1del <- -bit1.b1del 
bit3.b1del <- -bit1.b1del 
bit4.b1del <- bit1.b1del 

bit1.b2del <- c.copula2.be2del*d.n2
bit2.b2del <- -bit1.b2del 
bit3.b2del <- -bit1.b2del
bit4.b2del <- bit1.b2del  

bit2.thdel <- -bit1.thdel 
bit3.thdel <- -bit1.thdel 
bit4.thdel <- bit1.thdel 




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

  dl.ddelta <- weights*(y1.y2*c.copula.delta/p11+y1.cy2*(-c.copula.delta)/p10 + 
                       cy1.y2*(-c.copula.delta)/p01+cy1.cy2*c.copula.delta/p00)                          


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

  d2l.be1.del  <- -weights*(y1.y2*(bit1.b1del*p11-(c.copula.be1*d.n1*c.copula.delta))/p11^2+
                              y1.cy2*(bit2.b1del*p10-((1-c.copula.be1)*d.n1)*(-c.copula.delta))/p10^2+
                              cy1.y2*(bit3.b1del*p01-(-c.copula.be1*d.n1)*(-c.copula.delta))/p01^2+
                              cy1.cy2*(bit4.b1del*p00-((c.copula.be1-1)*d.n1)*c.copula.delta)/p00^2 )

  d2l.be2.del  <- -weights*(y1.y2*(bit1.b2del*p11-(c.copula.be2*d.n2*c.copula.delta))/p11^2+
                              y1.cy2*(bit2.b2del*p10-(-c.copula.be2*d.n2)*(-c.copula.delta))/p10^2+
                              cy1.y2*(bit3.b2del*p01-((1-c.copula.be2)*d.n2)*(-c.copula.delta))/p01^2+
                              cy1.cy2*(bit4.b2del*p00-((c.copula.be2-1)*d.n2)*c.copula.delta)/p00^2  )                           
                              
  d2l.del.del  <- -weights*(y1.y2*(bit1.del2*p11-c.copula.delta^2)/p11^2+
                              y1.cy2*(bit2.del2*p10-(-c.copula.delta)^2)/p10^2+
                              cy1.y2*(bit3.del2*p01-(-c.copula.delta)^2)/p01^2+
                              cy1.cy2*(bit4.del2*p00-c.copula.delta^2)/p00^2  ) 
                              
                              
 d2l.rho.del  <- -weights*(y1.y2*(bit1.thdel*p11-c.copula.theta*c.copula.delta)/p11^2+
                              y1.cy2*(bit2.thdel*p10-(c.copula.theta*c.copula.delta))/p10^2+
                              cy1.y2*(bit3.thdel*p01-(c.copula.theta*c.copula.delta))/p01^2+
                              cy1.cy2*(bit4.thdel*p00-c.copula.theta*c.copula.delta)/p00^2 ) 




  be1.be1 <- crossprod(X1*c(d2l.be1.be1),X1)
  be2.be2 <- crossprod(X2*c(d2l.be2.be2),X2)
  be1.be2 <- crossprod(X1*c(d2l.be1.be2),X2)
  be1.rho <- t(t(rowSums(t(X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(X2*c(d2l.be2.rho)))))
  be1.del <- t(t(rowSums(t(X1*c(d2l.be1.del)))))
  be2.del <- t(t(rowSums(t(X2*c(d2l.be2.del)))))


  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho, be1.del ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho, be2.del ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho), sum(d2l.rho.del) ),
              cbind( t(be1.del) , t(be2.del) , sum(d2l.rho.del), sum(d2l.del.del) )
            ) 
            
         res <- -sum(l.par)
         G   <- -c( colSums( c(dl.dbe1)*X1 ),
                    colSums( c(dl.dbe2)*X2 ),
                    sum( dl.drho ), sum( dl.ddelta )  )


if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0

     else{
        
    dimP1 <- dimP2 <- 0     
    S1 <- S2 <- matrix(0,1,1)   

    S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
    S <- do.call(adiag, lapply(S, unlist))

    ma1 <- matrix(0,gp1,gp1) 
    if(length(pPen1)!=0){ indP1 <- qu.mag$off[1]:(qu.mag$off[1]+qu.mag$rank[1]-1)
                          dimP1 <- length(indP1)
                          ma1[indP1,indP1] <- S[1:dimP1,1:dimP1]
                                } 
    ma2 <- matrix(0,gp2,gp2)
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
    
    if(lS1==1 && lS2==1) S.h <- adiag(ma1, ma2, 0, 0)
    if(lS1 >1 && lS2==1) S.h <- adiag(ma1, S1, ma2, 0, 0)
    if(lS1==1 && lS2 >1) S.h <- adiag(ma1, ma2, S2, 0, 0)
    if(lS1 >1 && lS2 >1) S.h <- adiag(ma1, S1, ma2, S2, 0, 0)
        
   
   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
   
         }


         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
         H   <- H + S.h  

         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, l.par=l.par, 
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho, dl.ddelta=dl.ddelta,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be1.del=d2l.be1.del,d2l.be2.del=d2l.be2.del,
              d2l.del.del=d2l.del.del,d2l.rho.del=d2l.rho.del ,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho,good=good, PL=PL, eqPL=eqPL,BivD=BivD)      

}




     
























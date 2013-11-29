bprobgHs.NPRE <- function(params, BivD, nC, nu, xi1, xi2, PL, eqPL, H.n=NULL, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, weights=weights, X1.d2, X2.d2, pPen1=NULL, pPen2=NULL, sp=NULL, qu.mag=NULL, gp1, gp2, fp, l.sp1, l.sp2, K, n, N, cuid, uidf, masses, NGQ=NULL, dat1all=NULL, dat2all=NULL, W=NULL){


  teta.st <- params[X1.d2+X2.d2+2*K+1]
  epsilon <- .Machine$double.eps*10^6
  
    if(BivD %in% c("N","T")      ){teta <- tanh(teta.st); if(teta %in% c(-1,1)) teta <- sign(teta)*0.9999999}
    if(BivD=="F")                  teta <- teta.st + epsilon
    if(BivD %in% c("C0", "C180") ) teta <- exp(teta.st) + epsilon
    if(BivD %in% c("C90","C270") ) teta <- -( exp(teta.st) + epsilon ) 
    if(BivD %in% c("J0", "J180") ) teta <- exp(teta.st) + 1 + epsilon 
    if(BivD %in% c("J90","J270") ) teta <- -( exp(teta.st) + 1 + epsilon ) 
    if(BivD %in% c("G0", "G180") ) teta <- exp(teta.st) + 1 
    if(BivD %in% c("G90","G270") ) teta <- -( exp(teta.st) + 1 )   


  l.par <- dl.dbe1 <- dl.dbe2 <- dl.drho <- d2l.be1.be1 <- d2l.be2.be2 <- d2l.be1.be2 <- d2l.be1.rho <-  d2l.be2.rho <- d2l.rho.rho <- matrix(0,nrow=n,ncol=K)
  
  
  for (u in 1:K){   

  eta1 <- X1(u)%*%params[1:(X1.d2+K)]
  eta2 <- X2(u)%*%params[(X1.d2+K+1):(X1.d2+X2.d2+2*K)]
  
  p1 <- pnorm(eta1); d.n1 <- dnorm(eta1)
  p2 <- pnorm(eta2); d.n2 <- dnorm(eta2) 

  criteria <- c(0,1)
  no.good <- apply(apply(cbind(p1,p2), c(1,2), `%in%`, criteria), 1, any)
  good <- no.good==FALSE

  p1 <- p1[good]
  p2 <- p2[good]
  d.n1 <- d.n1[good]
  d.n2 <- d.n2[good]
  eta1 <- eta1[good] 
  eta2 <- eta2[good] 
  y1.y2 <- y1.y2[good]
  y1.cy2 <- y1.cy2[good]
  cy1.y2 <- cy1.y2[good]
  cy1.cy2 <- cy1.cy2[good]
  weights <- weights[good]


  if(BivD=="N") C.copula <- pmax( abs(pbinorm( qnorm(p1), qnorm(p2), cov12=teta)), 1000*.Machine$double.eps ) else C.copula <- BiCopCDF(p1,p2, nC, par=teta, par2=nu)

  
  p11 <- pmax( C.copula, 1000*.Machine$double.eps )
  p10 <- pmax( p1 - p11, 1000*.Machine$double.eps )
  p01 <- pmax( p2 - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )
  
  
  l.par[,u] <- weights*( y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) )
  
  dH <- copgHs(p1,p2,eta1=NULL,eta2=NULL,teta,teta.st,xi1=NULL,xi1.st=NULL,xi2=NULL,xi2.st=NULL,BivD,nC,nu,PL,eqPL)
  
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
  
  bit1.th2 <- dH$bit1.th2
  bit2.th2 <- -bit1.th2 
  bit3.th2 <- -bit1.th2 
  bit4.th2 <- bit1.th2 
  
  
  dl.dbe1[,u] <-  weights*d.n1*( (y1.y2*c.copula.be1/p11)  +
                      (y1.cy2*(1-c.copula.be1)/p10) +
                      (cy1.y2*c.copula.be1/(-p01)) +
                      (cy1.cy2*(c.copula.be1-1)/p00) )
                                
  dl.dbe2[,u] <-  weights*d.n2*( (y1.y2*c.copula.be2/p11)  +
                            (y1.cy2*c.copula.be2/(-p10)) +
                                (cy1.y2*(1-c.copula.be2)/(p01)) +
                                (cy1.cy2*(c.copula.be2-1)/p00) )

  dl.drho[,u] <- weights*( y1.y2*c.copula.theta/p11+y1.cy2*(-c.copula.theta)/p10 + 
                       cy1.y2*(-c.copula.theta)/p01+cy1.cy2*c.copula.theta/p00 ) 


  d2l.be1.be1[,u]  <- -weights*(y1.y2*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2+
                              y1.cy2*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2+
                              cy1.y2*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2+
                              cy1.cy2*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2 )

  d2l.be2.be2[,u]  <- -weights*(y1.y2*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              y1.cy2*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2+
                              cy1.y2*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2+
                              cy1.cy2*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2 )

  d2l.be1.be2[,u]  <- -weights*(y1.y2*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2+
                              y1.cy2*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2+
                              cy1.y2*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2+
                              cy1.cy2*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2 )

  d2l.be1.rho[,u]  <- -weights*(y1.y2*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )

  d2l.be2.rho[,u]  <- -weights*(y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )

  d2l.rho.rho[,u]  <- -weights*(y1.y2*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                              y1.cy2*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2+
                              cy1.y2*(bit3.th2*p01-(-c.copula.theta)^2)/p01^2+
                              cy1.cy2*(bit4.th2*p00-c.copula.theta^2)/p00^2 ) 

  } 

  Wp1 <- exp(l.par)
  Wp2 <- matrix(0,nrow=N,ncol=K)
  for (i in 1:N) if (uidf[i]>1) {Wp2[i,] <- apply(Wp1[(cuid[i]+1):(cuid[i+1]),],2,prod)
                                 }else{Wp2[i,] <- Wp1[(cuid[i]+1):(cuid[i+1]),]}
  Wp3 <- t(masses*t(Wp2))
  W   <- Wp3/apply(Wp3,1,sum)
  We  <- matrix(rep(c(W),rep(uidf,K)),ncol=K)

  be1.be1 <- be2.be2 <- be1.be2 <- be1.rho <- be2.rho <- rho.rho <- g1 <- g2 <- g3 <- 0 
  
  for (u in 1:K){  
  
  X1u <- X1(u)[good,]
  X2u <- X2(u)[good,]

  g1 <- g1 + colSums( c(dl.dbe1[,u])*We[,u]*X1u )
  g2 <- g2 + colSums( c(dl.dbe2[,u])*We[,u]*X2u )
  g3 <- g3 + sum( c(dl.drho[,u])*We[,u] )
  be1.be1 <- be1.be1 + crossprod(X1u*c(d2l.be1.be1[,u])*We[,u],X1u)
  be2.be2 <- be2.be2 + crossprod(X2u*c(d2l.be2.be2[,u])*We[,u],X2u)
  be1.be2 <- be1.be2 + crossprod(X1u*c(d2l.be1.be2[,u])*We[,u],X2u)
  be1.rho <- be1.rho + t(t(rowSums(t(X1u*c(d2l.be1.rho[,u])*We[,u]))))
  be2.rho <- be2.rho + t(t(rowSums(t(X2u*c(d2l.be2.rho[,u])*We[,u]))))
  rho.rho <- rho.rho + c(d2l.rho.rho[,u])*We[,u] 
  }
    
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ),
              cbind( t(be1.be2) , be2.be2    , be2.rho ),
              cbind( t(be1.rho) , t(be2.rho) , sum(rho.rho) )
            )

  G      <- -c(g1,g2,g3)
  res    <- -sum(We*l.par)
  masses <- apply(W,2,sum)/sum(W)

  #################################################### 
  
#if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0
#
#     else{
#
#       S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
#       S <- do.call(adiag, lapply(S, unlist)) 
#
#    if(l.sp1!=0 && l.sp2!=0) S.h <- adiag(matrix(0,K+gp1-1,K+gp1-1),
#                                          S[1:(X1.d2-gp1),1:(X1.d2-gp1)],
#                                          matrix(0,K+gp2-1,K+gp2-1),
#                                          S[(X1.d2-(gp1-1)):dim(S)[2],(X1.d2-(gp1-1)):dim(S)[2]],
#                                          0)
#
#    if(l.sp1==0 && l.sp2!=0) S.h <- adiag(matrix(0,K+gp1-1,K+gp1-1), matrix(0,K+gp2-1,K+gp2-1), S, 0)
#    if(l.sp1!=0 && l.sp2==0) S.h <- adiag(matrix(0,K+gp1-1,K+gp1-1), S, matrix(0,K+gp2-1,K+gp2-1), 0)
#
#                      
#   S.h1 <- 0.5*crossprod(params,S.h)%*%params
#   S.h2 <- S.h%*%params
#         }
         
         
                   
if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0

     else{
        
    dimP1 <- dimP2 <- 0     
    S1 <- S2 <- matrix(0,1,1)   

    S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
    S <- do.call(adiag, lapply(S, unlist))

    ma1 <- matrix(0,gp1+(K-1),gp1+(K-1)) 
    if(length(pPen1)!=0){ indP1 <- qu.mag$off[1]:(qu.mag$off[1]+qu.mag$rank[1]-1)
                          dimP1 <- length(indP1)
                          ma1[indP1,indP1] <- S[1:dimP1,1:dimP1]
                                } 
    ma2 <- matrix(0,gp2+(K-1),gp2+(K-1))
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
    
    if(lS1==1 && lS2==1) S.h <- adiag(ma1, ma2, 0)
    if(lS1 >1 && lS2==1) S.h <- adiag(ma1, S1, ma2, 0)
    if(lS1==1 && lS2 >1) S.h <- adiag(ma1, ma2, S2, 0)
    if(lS1 >1 && lS2 >1) S.h <- adiag(ma1, S1, ma2, S2, 0)
        
   
   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
   
         }         
         
 

         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
         H   <- H + S.h  

         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, masses=masses,
              eta1=eta1,#X1[,-c(1:K)]%*%params[(K+1):(X1.d2+K)], 
              eta2=eta2,#X2[,-c(1:K)]%*%params[(X1.d2+2*K+1):(X1.d2+X2.d2+2*K)],
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho, We=We, good=good, PL=PL,eqPL=eqPL,BivD=BivD)



}

























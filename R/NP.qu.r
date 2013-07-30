NP.qu <- function(SemiParFit, y1, y2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, X1.d2, X2.d2, qu.mag=NULL, gp1, gp2, 
                    fp, l.sp1, l.sp2, weights, K, n, N, cuid, uidf, NGQ=NULL, dat1all=NULL, dat2all=NULL, W=NULL){
                            
        #T.sv <- bprobNP.H(SemiParFit$fit$argument, y1=y1, y2=y2, 
        #                 y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1, 
        #                 X1=X1, X2=X2, 
        #                 X1.d2=X1.d2, X2.d2=X2.d2, sp=SemiParFit$sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, weights=weights,
        #                 K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=SemiParFit$masses, NGQ=NGQ, dat1all=dat1all, dat2all=dat2all, W=W)
      
      
  params <- SemiParFit$fit$argument   
  masses <- SemiParFit$fit$masses
  sp <- SemiParFit$sp
      
  corr.st <- params[(X1.d2+X2.d2+1)]
  corr    <- tanh(corr.st)
  d.r    <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr^2) )  
  drh.drh.st <- 1/cosh(corr.st)^2

  l.par <- dl.dbe1 <- dl.dbe2 <- dl.drho <- d2l.be1.be1 <- d2l.be2.be2 <- d2l.be1.be2 <- d2l.be1.rho <-  d2l.be2.rho <- d2l.rho.rho <- matrix(0,nrow=n,ncol=K)
  
  for (u in 1:K){   

  eta1 <- X1(u)%*%params[1:(X1.d2+K)]
  eta2 <- X2(u)%*%params[(X1.d2+K+1):(X1.d2+X2.d2+2*K)]

  p1 <- pnorm(eta1)
  p2 <- pnorm(eta2)
  
  C.copula <- pnorm2( eta1, eta2, cov12=corr)
  
  p11 <- pmax( C.copula, 1000*.Machine$double.eps )
  p10 <- pmax( p1 - p11, 1000*.Machine$double.eps )
  p01 <- pmax( p2 - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )
  
  d.n1   <- dnorm(eta1) 
  d.n2   <- dnorm(eta2) 
  d.n1n2 <- dnorm2(eta1,eta2, rho=corr) 

  A   <- pnorm( (eta2-corr*eta1)*d.r )
  B   <- pnorm( (eta1-corr*eta2)*d.r )

  A   <- pnorm( (eta2-corr*eta1)*d.r )
  A.c <- 1 - A
  B   <- pnorm( (eta1-corr*eta2)*d.r )
  B.c <- 1 - B

  l.par[,u] <- weights*( y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) )
  
  
    c.copula.be1 <- A  
    c.copula.be2 <- B 
    c.copula.theta <-  d.n1n2*drh.drh.st 
    
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
                                                
  
  
 
c.copula2.be1 <- dnorm((eta2-corr*eta1)*d.r)*-corr*d.r*sqrt(2*pi)/exp(-eta1^2/2)  
c.copula2.be2 <- dnorm((eta1-corr*eta2)*d.r)*-corr*d.r*sqrt(2*pi)/exp(-eta2^2/2)
             
bit1.b1b1 <- c.copula2.be1*(d.n1)^2-c.copula.be1*d.n1*eta1
bit2.b1b1 <- -d.n1*eta1-bit1.b1b1
bit3.b1b1 <- -bit1.b1b1
bit4.b1b1 <- -bit2.b1b1

bit1.b2b2 <- c.copula2.be2*(d.n2)^2-c.copula.be2*d.n2*eta2
bit2.b2b2 <- -bit1.b2b2
bit3.b2b2 <- -d.n2*eta2-bit1.b2b2
bit4.b2b2 <- -bit3.b2b2
##

##
c.copula2.be1be2 <- d.r*exp( -(  corr^2*(eta1^2+eta2^2)-2*corr*eta1*eta2     )/(2*(1-corr^2))          )

bit1.b1b2 <- c.copula2.be1be2 * d.n1 *d.n2
bit2.b1b2 <- -bit1.b1b2
bit3.b1b2 <- -bit1.b1b2
bit4.b1b2 <- bit1.b1b2
##

##
c.copula2.be1th <- -(dnorm((eta2 - tanh(corr.st) * eta1)/sqrt(1 - tanh(corr.st)^2)) * 
    (1/cosh(corr.st)^2 * eta1/sqrt(1 - tanh(corr.st)^2) - (eta2 - 
        tanh(corr.st) * eta1) * (0.5 * (2 * (1/cosh(corr.st)^2 * 
        tanh(corr.st)) * (1 - tanh(corr.st)^2)^-0.5))/sqrt(1 - 
        tanh(corr.st)^2)^2)) 

        
bit1.b1th <- c.copula2.be1th*d.n1

bit2.b1th <- -bit1.b1th 
bit3.b1th <- -bit1.b1th 
bit4.b1th <- bit1.b1th 
##

##
c.copula2.be2th <- -(dnorm((eta1 - tanh(corr.st) * eta2)/sqrt(1 - tanh(corr.st)^2)) * 
    (1/cosh(corr.st)^2 * eta2/sqrt(1 - tanh(corr.st)^2) - (eta1 - 
        tanh(corr.st) * eta2) * (0.5 * (2 * (1/cosh(corr.st)^2 * 
        tanh(corr.st)) * (1 - tanh(corr.st)^2)^-0.5))/sqrt(1 - 
        tanh(corr.st)^2)^2))


bit1.b2th <- c.copula2.be2th*d.n2
bit2.b2th <- -bit1.b2th 
bit3.b2th <- -bit1.b2th 
bit4.b2th <- bit1.b2th 
##

##
bit1.th2 <- (2 * pi * (0.5 * (2 * (1/cosh(corr.st)^2 * tanh(corr.st)) * (1 - 
    tanh(corr.st)^2)^-0.5))/(2 * pi * sqrt(1 - tanh(corr.st)^2))^2 * 
    exp(-1/(2 * (1 - tanh(corr.st)^2)) * (eta1^2 + eta2^2 - 2 * 
        tanh(corr.st) * eta1 * eta2)) - 1/(2 * pi * sqrt(1 - 
    tanh(corr.st)^2)) * (exp(-1/(2 * (1 - tanh(corr.st)^2)) * 
    (eta1^2 + eta2^2 - 2 * tanh(corr.st) * eta1 * eta2)) * (-1/(2 * 
    (1 - tanh(corr.st)^2)) * (2 * (1/cosh(corr.st)^2) * eta1 * 
    eta2) + 2 * (2 * (1/cosh(corr.st)^2 * tanh(corr.st)))/(2 * 
    (1 - tanh(corr.st)^2))^2 * (eta1^2 + eta2^2 - 2 * tanh(corr.st) * 
    eta1 * eta2))))/cosh(corr.st)^2 - 1/(2 * pi * sqrt(1 - tanh(corr.st)^2)) * 
    exp(-1/(2 * (1 - tanh(corr.st)^2)) * (eta1^2 + eta2^2 - 2 * 
        tanh(corr.st) * eta1 * eta2)) * 1 * (2 * (sinh(corr.st) * 
    cosh(corr.st)))/(cosh(corr.st)^2)^2 

bit2.th2 <- -bit1.th2 
bit3.th2 <- -bit1.th2 
bit4.th2 <- bit1.th2 







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
  g1 <- g1 - colSums( c(dl.dbe1[,u])*We[,u]*X1(u) )
  g2 <- g2 - colSums( c(dl.dbe2[,u])*We[,u]*X2(u) )
  g3 <- g3 - sum( c(dl.drho[,u])*We[,u] )
  be1.be1 <- be1.be1 + crossprod(X1(u)*c(d2l.be1.be1[,u])*We[,u],X1(u))
  be2.be2 <- be2.be2 + crossprod(X2(u)*c(d2l.be2.be2[,u])*We[,u],X2(u))
  be1.be2 <- be1.be2 + crossprod(X1(u)*c(d2l.be1.be2[,u])*We[,u],X2(u))
  be1.rho <- be1.rho + t(t(rowSums(t(X1(u)*c(d2l.be1.rho[,u])*We[,u]))))
  be2.rho <- be2.rho + t(t(rowSums(t(X2(u)*c(d2l.be2.rho[,u])*We[,u]))))
  rho.rho <- rho.rho + c(d2l.rho.rho[,u])*We[,u] 
  }
    
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ),
              cbind( t(be1.be2) , be2.be2    , be2.rho ),
              cbind( t(be1.rho) , t(be2.rho) , sum(rho.rho) )
            )

  G      <- c(g1,g2,g3)
  res    <- -sum(We*l.par)
  masses <- apply(W,2,sum)/sum(W)

  #################################################### 
  
if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0

     else{

       S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
       S <- do.call(adiag, lapply(S, unlist)) 

    if(l.sp1!=0 && l.sp2!=0) S.h <- adiag(matrix(0,K+gp1-1,K+gp1-1),
                                          S[1:(X1.d2-gp1),1:(X1.d2-gp1)],
                                          matrix(0,K+gp2-1,K+gp2-1),
                                          S[(X1.d2-(gp1-1)):dim(S)[2],(X1.d2-(gp1-1)):dim(S)[2]],
                                          0)

    if(l.sp1==0 && l.sp2!=0) S.h <- adiag(matrix(0,K+gp1-1,K+gp1-1), matrix(0,K+gp2-1,K+gp2-1), S, 0)
    if(l.sp1!=0 && l.sp2==0) S.h <- adiag(matrix(0,K+gp1-1,K+gp1-1), S, matrix(0,K+gp2-1,K+gp2-1), 0)

                      
   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
         }

         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
         H   <- H + S.h  

      T.sv <-   list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, masses=masses,
              eta1=eta1,  
              eta2=eta2,  
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho, We=We, W=W,l.par=l.par,Wp3=Wp3)


       
       
       
       
        #T.sv <- SemiParFit$fit                 

        npar <- K + X1.d2 + K + X2.d2 + 1 + K - 1
        H.cor <- G.cor <- matrix(0,nrow=npar,ncol=npar)

        for (i in 1:N){
              h <- matrix(0,nrow=npar,ncol=1)
           for (l in 1:K){

              c1 <- T.sv$dl.dbe1[(cuid[i]+1):(cuid[i+1]),l]*X1(l)[(cuid[i]+1):(cuid[i+1]),]
              if (uidf[i]>1) c1<-apply(c1,2,sum)

              c2 <- T.sv$dl.dbe2[(cuid[i]+1):(cuid[i+1]),l]*X2(l)[(cuid[i]+1):(cuid[i+1]),]
              if (uidf[i]>1) c2<-apply(c2,2,sum)

              c3 <- sum(T.sv$dl.drho[(cuid[i]+1):(cuid[i+1]),l])

              if (l < K) {c4 <- array(0,K-1); c4[l] <- (1/T.sv$masses[l])}else{ c4 <- array((-1/T.sv$masses[K]),K-1)}

               g <- as.matrix(c(c1,c2,c3,c4))

              H.cor <- H.cor + T.sv$W[i,l]*tcrossprod(g)

              h <- h + T.sv$W[i,l]*g 
           }
           G.cor <- G.cor + tcrossprod(h)
        }

        I.prob1 <- apply(t((1/T.sv$masses^2)*t(T.sv$W)),2,sum)
        if(K!=2) I.prob  <- diag(I.prob1[1:K-1]) + matrix(I.prob1[K],K-1,K-1) else I.prob  <- sum(I.prob1) # I.prob1[1:K-1] + I.prob1[K] 

        He <- adiag(T.sv$hessian,I.prob) - H.cor + G.cor   

        logLik <- sum(log(apply(T.sv$Wp3,1,sum))) 

        u1 <- SemiParFit$fit$argument[1:K]
        u2 <- SemiParFit$fit$argument[(K+X1.d2+1):(K+X1.d2+K)]

        nw <- T.sv$Wp3/matrix(rep(apply(T.sv$Wp3,1,sum),K),ncol=K)

        eb.u1 <- apply(t(u1*t(nw)),1,sum)
        eb.u2 <- apply(t(u2*t(nw)),1,sum)

        Eb.u1 <- rep(eb.u1,uidf)
        Eb.u2 <- rep(eb.u2,uidf)

    list(He=He,logLik=logLik,eb.u1=eb.u1,eb.u2=eb.u2,Eb.u1=Eb.u1,Eb.u2=Eb.u2,T.sv=T.sv)

    }


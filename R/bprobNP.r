bprobNP <- function(params, BivD=NULL, nC=NULL, nu=NULL, H.n=NULL, y1, y2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, weights=weights, X1.d2, X2.d2, pPen1=NULL, pPen2=NULL, sp=NULL, qu.mag=NULL, gp1, gp2, fp, l.sp1, l.sp2, K, n, N, cuid, uidf, masses, NGQ=NULL, dat1all=NULL, dat2all=NULL, W=NULL){

  corr.st <- params[X1.d2+X2.d2+2*K+1]
  corr    <- tanh(corr.st)
  d.r     <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr^2) )
  drh.drh.st   <- 4*exp(2*corr.st)/(exp(2*corr.st)+1)^2

  l.par <- dl.dbe1 <- dl.dbe2 <- dl.drho <- d2l.be1.be1 <- d2l.be2.be2 <- d2l.be1.be2 <- d2l.be1.rho <-  d2l.be2.rho <- d2l.rho.rho <- matrix(0,nrow=n,ncol=K)
  
  for (u in 1:K){   

  eta1 <- X1(u)%*%params[1:(X1.d2+K)]
  eta2 <- X2(u)%*%params[(X1.d2+K+1):(X1.d2+X2.d2+2*K)]

  A   <- pnorm( (eta2-corr*eta1)*d.r )
  A.c <- 1 - A
  B   <- pnorm( (eta1-corr*eta2)*d.r )
  B.c <- 1 - B

  p11 <- pmax( abs(pnorm2( eta1, eta2, cov12=corr)), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
  p01 <- pmax( pnorm(eta2) - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )

  d.n1   <- dnorm(eta1)
  d.n2   <- dnorm(eta2)
  d.n1n2 <- dnorm2(eta1,eta2,rho=corr)

  l.par[,u] <- weights*( y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) )
  dl.dbe1[,u] <- weights*d.n1*((y1.y2/p11-cy1.y2/p01)*A+(y1.cy2/p10-cy1.cy2/p00)*A.c)
  dl.dbe2[,u] <- weights*d.n2*((y1.y2/p11-y1.cy2/p10)*B+(cy1.y2/p01-cy1.cy2/p00)*B.c)
  dl.drho[,u] <- weights*d.n1n2*(y1.y2/p11-y1.cy2/p10-cy1.y2/p01+cy1.cy2/p00)*drh.drh.st
  d2l.be1.be1[,u]  <- -weights*(d.n1^2*(A^2*(-1/p11-1/p01)+A.c^2*(-1/p10-1/p00)))
  d2l.be2.be2[,u]  <- -weights*(d.n2^2*(B^2*(-1/p11-1/p10)+B.c^2*(-1/p01-1/p00)))
  d2l.be1.be2[,u]  <- -weights*(d.n1*d.n2*(A*B.c/p01-A*B/p11+A.c*B/p10-A.c*B.c/p00))
  d2l.be1.rho[,u]  <- -weights*(-d.n1*d.n1n2*(A*(1/p11+1/p01)-A.c*(1/p10+1/p00)))*drh.drh.st
  d2l.be2.rho[,u]  <- -weights*(-d.n2*d.n1n2*(B*(1/p11+1/p10)-B.c*(1/p01+1/p00)))*drh.drh.st
  d2l.rho.rho[,u]  <- -weights*(-d.n1n2^2*(1/p11+1/p01+1/p10+1/p00))*drh.drh.st^2

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

         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, masses=masses,
              eta1=eta1,#X1[,-c(1:K)]%*%params[(K+1):(X1.d2+K)], 
              eta2=eta2,#X2[,-c(1:K)]%*%params[(X1.d2+2*K+1):(X1.d2+X2.d2+2*K)],
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho, We=We)



}

























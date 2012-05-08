bprobNP.H <- function(params, dat, dat1, dat2, dat1p=NULL, dat2p=NULL, X1.d2, X2.d2, S=NULL, gam1, gam2, fp, K, n, N, cuid, uidf, masses){

  q1 <- 2*dat[,1]-1
  q2 <- 2*dat[,2]-1
  
  corr    <- tanh(params[X1.d2+X2.d2+2*K+1])
  corr.sq <- q1*q2*corr
  delta <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr.sq ^2) )

  l.par <- dl.dbe1 <- dl.dbe2 <- dl.drho <- d2l.be1.be1 <- d2l.be2.be2 <- d2l.be1.be2 <- d2l.be1.rho <-  d2l.be2.rho <- d2l.rho.rho <- matrix(0,nrow=n,ncol=K)
  
  for (u in 1:K){   

  eta1 <- dat1(u)%*%params[1:(X1.d2+K)]
  eta2 <- dat2(u)%*%params[(X1.d2+K+1):(X1.d2+X2.d2+2*K)]
  w1 <- q1*eta1
  w2 <- q2*eta2
  
  g1 <- dnorm(w1)*pnorm( (w2-corr.sq*w1) * delta )
  g2 <- dnorm(w2)*pnorm( (w1-corr.sq*w2) * delta )
  
  PHI2   <- pmax(pnorm2( w1, w2, corr.sq),1000*.Machine$double.eps)
  L.PHI2 <- log(PHI2) 
  phi2   <- dnorm2( w1, w2, corr.sq)
       
  v1 <- delta*( w2 - corr.sq*w1 )
  v2 <- delta*( w1 - corr.sq*w2 )

  l.par[,u] <- L.PHI2
  dl.dbe1[,u] <- q1*g1/PHI2
  dl.dbe2[,u] <- q2*g2/PHI2
  dl.drho[,u] <- q1*q2*phi2/PHI2
  d2l.be1.be1[,u]  <- -(-w1*g1 - corr.sq*phi2 - g1^2/PHI2)/PHI2
  d2l.be2.be2[,u]  <- -(-w2*g2 - corr.sq*phi2 - g2^2/PHI2)/PHI2
  d2l.be1.be2[,u]  <- -(phi2/PHI2 - g1*g2/PHI2^2)*q1*q2
  d2l.be1.rho[,u]  <- -(corr.sq*delta*v1 - w1 - g1/PHI2)*phi2/PHI2*q2
  d2l.be2.rho[,u]  <- -(corr.sq*delta*v2 - w2 - g2/PHI2)*phi2/PHI2*q1
  d2l.rho.rho[,u]  <- -( phi2/PHI2*( delta^2*corr.sq*( 1 - delta^2*(w1^2+w2^2 - 2*corr.sq*w1*w2) ) + delta^2*w1*w2 - phi2/PHI2 ) )

  } 

  Wp1 <- exp(l.par)
  Wp2 <- matrix(0,nrow=N,ncol=K)
  for (i in 1:N) if (uidf[i]>1) {Wp2[i,] <- apply(Wp1[(cuid[i]+1):(cuid[i+1]),],2,prod)
                                 }else{Wp2[i,] <- Wp1[(cuid[i]+1):(cuid[i+1]),]}
  Wp3 <- t(masses*t(Wp2))

  W <- Wp3/apply(Wp3,1,sum)
  We <- matrix(rep(c(W),rep(uidf,K)),ncol=K)

  be1.be1 <- be2.be2 <- be1.be2 <- be1.rho <- be2.rho <- rho.rho <- g1 <- g2 <- g3 <- 0 
  
  for (u in 1:K){  
  g1 <- g1 - colSums( c(dl.dbe1[,u])*We[,u]*dat1(u) )
  g2 <- g2 - colSums( c(dl.dbe2[,u])*We[,u]*dat2(u) )
  g3 <- g3 - sum( c(dl.drho[,u])*We[,u] )
  be1.be1 <- be1.be1 + crossprod(dat1(u)*c(d2l.be1.be1[,u])*We[,u],dat1(u))
  be2.be2 <- be2.be2 + crossprod(dat2(u)*c(d2l.be2.be2[,u])*We[,u],dat2(u))
  be1.be2 <- be1.be2 + crossprod(dat1(u)*c(d2l.be1.be2[,u])*We[,u],dat2(u))
  be1.rho <- be1.rho + t(t(rowSums(t(dat1(u)*c(d2l.be1.rho[,u])*We[,u]))))
  be2.rho <- be2.rho + t(t(rowSums(t(dat2(u)*c(d2l.be2.rho[,u])*We[,u]))))
  rho.rho <- rho.rho + c(d2l.rho.rho[,u])*We[,u] 
  }
    
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ),
              cbind( t(be1.be2) , be2.be2    , be2.rho ),
              cbind( t(be1.rho) , t(be2.rho) , sum(rho.rho) )
            )

  G      <- c(g1,g2,g3)
  res    <- -sum(We*l.par)
  masses <- apply(W,2,sum)/sum(W)
  
  
   if( ( length(gam1$smooth)==0 && length(gam2$smooth)==0 ) || fp==TRUE){
  
         list(value=res, gradient=G, hessian=H, l=res, masses=masses,
              dl.dbe1=dl.dbe1,dl.dbe2=dl.dbe2,dl.drho=dl.drho,W=W,l.par=l.par,Wp3=Wp3)
         
  }else{
         S.h <- adiag(matrix(0,K+gam1$nsdf-1,K+gam1$nsdf-1),
                      S[1:(X1.d2-gam1$nsdf),1:(X1.d2-gam1$nsdf)],
                      matrix(0,K+gam2$nsdf-1,K+gam2$nsdf-1),
                      S[(X1.d2-(gam1$nsdf-1)):dim(S)[2],(X1.d2-(gam1$nsdf-1)):dim(S)[2]],
                      0)
          
         S.res <- res 
         res <- S.res + 0.5*crossprod(params,S.h)%*%params
         G   <- G + S.h%*%params
         H   <- H + S.h
         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, masses=masses,
              dl.dbe1=dl.dbe1,dl.dbe2=dl.dbe2,dl.drho=dl.drho,W=W,l.par=l.par,Wp3=Wp3)
   }
}


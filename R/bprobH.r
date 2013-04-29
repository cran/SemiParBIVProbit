bprobH <- function(params, y1, y2, q1, q2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, weights=weights, X1.d2, X2.d2, sp=NULL, qu.mag=NULL, gp1, gp2, fp=NULL, l.sp1=NULL, l.sp2=NULL, K=NULL, n=NULL, N=NULL, cuid=NULL, uidf=NULL, masses=NULL, NGQ=NULL, dat1all=NULL, dat2all=NULL, W=NULL){

  eta1 <- X1%*%params[1:X1.d2]
  eta2 <- X2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  corr.st <- params[(X1.d2+X2.d2+1)]
  corr    <- tanh(corr.st)
  drh.drh.st  <- 4*exp(2*corr.st)/(exp(2*corr.st)+1)^2
  drh.drh.st2 <- 8*exp(2*corr.st)*(1-exp(2*corr.st)) / (exp(2*corr.st)+1)^3
  corr.sq <- q1*q2*corr
  delta <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr.sq^2) )
  
  w1 <- q1*eta1
  w2 <- q2*eta2
  v1 <- delta*( w2 - corr.sq*w1 )
  v2 <- delta*( w1 - corr.sq*w2 )  
  
  g1 <- dnorm(w1)*pnorm( v1 )
  g2 <- dnorm(w2)*pnorm( v2 )
  
  PHI2   <- pmax(abs(pnorm2( w1, w2, cov12=corr.sq)),1000*.Machine$double.eps)
  phi2   <- dnorm2( w1, w2, rho=corr.sq)

  l.par <- weights*log(PHI2) 
  
  dl.dbe1 <- weights*q1*g1/PHI2
  dl.dbe2 <- weights*q2*g2/PHI2
  dl.drho <- weights*q1*q2*phi2/PHI2*drh.drh.st
  
  d2l.be1.be1 <- -weights*(-w1*g1 - corr.sq*phi2 - g1^2/PHI2)/PHI2  
  d2l.be2.be2 <- -weights*(-w2*g2 - corr.sq*phi2 - g2^2/PHI2)/PHI2 
  d2l.be1.be2 <- -weights*(q1*q2)/PHI2*(phi2 - g1*g2/PHI2) 
  d2l.be1.rho <- -weights*(corr.sq*delta*v1 - w1 - g1/PHI2)*q2*phi2/PHI2*drh.drh.st
  d2l.be2.rho <- -weights*(corr.sq*delta*v2 - w2 - g2/PHI2)*q1*phi2/PHI2*drh.drh.st
  d2l.rho.rho <- -weights*( phi2/PHI2*( delta^2*corr.sq*( 1 - delta^2*(w1^2+w2^2 - 2*corr.sq*w1*w2) ) + delta^2*w1*w2 - phi2/PHI2 )*drh.drh.st^2 + q1*q2*phi2/PHI2*drh.drh.st2 )
 
  be1.be1 <- crossprod(X1*c(d2l.be1.be1),X1)
  be2.be2 <- crossprod(X2*c(d2l.be2.be2),X2)
  be1.be2 <- crossprod(X1*c(d2l.be1.be2),X2)
  be1.rho <- t(t(rowSums(t(X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(X2*c(d2l.be2.rho)))))
  
  p11 <- pmax( pnorm2( eta1, eta2, cov12=corr), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
  p01 <- pmax( pnorm(eta2) - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )

  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
            ) 

         res <- -sum(l.par)
         G   <- c( -colSums( c(dl.dbe1)*X1 ),
                   -colSums( c(dl.dbe2)*X2 ),
                   -sum( dl.drho )  )

         list(value=res, gradient=G, hessian=H, l=res, 
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho) 
      

}


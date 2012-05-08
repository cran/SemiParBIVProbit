bprobSSH <- function(params, dat, dat1, dat2, dat1p=NULL, dat2p=NULL, X1.d2, X2.d2, S=NULL, gam1, gam2, fp, K=NULL, n=NULL, N=NULL, cuid=NULL, uidf=NULL, masses=NULL){

    q1 <- 2*dat[,1]-1
    q2 <- 2*dat[,2]-1 
    eta1 <- dat1%*%params[1:X1.d2]
    eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
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
    
    PHI2   <- pmax(pnorm2( w1, w2, corr.sq),1000*.Machine$double.eps)
    phi2   <- dnorm2( w1, w2, corr.sq)
  
  
    y1  <- dat[,1]
    cy1 <- (1-y1)  	
    d.n1 <- dnorm(eta1) 
    p0   <- pmax( pnorm(-eta1), 1000*.Machine$double.eps )
    CY1  <- cy1*(-(d.n1/p0)^2 + d.n1*eta1/p0) 
  
  
    l.par <- y1*log(PHI2) + cy1*log(p0)  
    
    dl.dbe1 <- y1*q1*g1/PHI2 - d.n1*cy1/p0
    dl.dbe2 <- y1*q2*g2/PHI2
    dl.drho <- y1*q1*q2*phi2/PHI2*drh.drh.st
    
    d2l.be1.be1 <- -y1*(-w1*g1 - corr.sq*phi2 - g1^2/PHI2)/PHI2 - CY1  
    d2l.be2.be2 <- -y1*(-w2*g2 - corr.sq*phi2 - g2^2/PHI2)/PHI2 
    d2l.be1.be2 <- -y1*(q1*q2)/PHI2*(phi2 - g1*g2/PHI2) 
    d2l.be1.rho <- -y1*(corr.sq*delta*v1 - w1 - g1/PHI2)*q2*phi2/PHI2*drh.drh.st
    d2l.be2.rho <- -y1*(corr.sq*delta*v2 - w2 - g2/PHI2)*q1*phi2/PHI2*drh.drh.st
    d2l.rho.rho <- -y1*( phi2/PHI2*( delta^2*corr.sq*( 1 - delta^2*(w1^2+w2^2 - 2*corr.sq*w1*w2) ) + delta^2*w1*w2 - phi2/PHI2 )*drh.drh.st^2 + q1*q2*phi2/PHI2*drh.drh.st2 )
 
             
  p11 <- pmax( pnorm2( eta1, eta2, corr), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
           
       
  be1.be1 <- crossprod(dat1*c(d2l.be1.be1),dat1)
  be2.be2 <- crossprod(dat2*c(d2l.be2.be2),dat2)
  be1.be2 <- crossprod(dat1*c(d2l.be1.be2),dat2)
  be1.rho <- t(t(rowSums(t(dat1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(dat2*c(d2l.be2.rho)))))
  
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
            ) 
            
         res <- -sum(l.par)
         G   <- c( -colSums( c(dl.dbe1)*dat1 ),
                   -colSums( c(dl.dbe2)*dat2 ),
                   -sum( dl.drho )  )            

         list(value=res, gradient=G, hessian=H, l=res, 
              p11=p11, p10=p10, p0=p0, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho) 
     

}


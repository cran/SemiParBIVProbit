bprobSS <- function(params, dat, dat1, dat2, dat1p=NULL, dat2p=NULL, X1.d2, X2.d2, S=NULL, gam1, gam2, fp, K=NULL, n=NULL, N=NULL, cuid=NULL, uidf=NULL, masses=NULL){

  eta1 <- dat1%*%params[1:X1.d2]
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  corr.st <- params[(X1.d2+X2.d2+1)]
  corr    <- tanh(corr.st)

  y1.y2  <- dat[,1]*dat[,2]
  y1.cy2 <- dat[,1]*(1-dat[,2])
  cy1    <- (1-dat[,1])

  d.r  <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr^2) )

  A   <- pnorm( (eta2-corr*eta1)*d.r )
  A.c <- 1 - A
  B   <- pnorm( (eta1-corr*eta2)*d.r )

  p11 <- pmax( pnorm2( eta1, eta2, corr), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
  p0  <- pmax( pnorm(-eta1), 1000*.Machine$double.eps )

  d.n1   <- dnorm(eta1) 
  d.n2   <- dnorm(eta2) 
  d.n1n2 <- dnorm2(eta1,eta2,corr) 

  l.par <- y1.y2*log(p11)+y1.cy2*log(p10)+cy1*log(p0) 

  drh.drh.st   <- 4*exp(2*corr.st)/(exp(2*corr.st)+1)^2
  
  dl.dbe1 <- d.n1*( y1.y2/p11*A  + y1.cy2/p10*A.c - cy1/p0 )  
  dl.dbe2 <- d.n2*B*( y1.y2/p11 - y1.cy2/p10)  
  dl.drho <- d.n1n2*(y1.y2/p11 - y1.cy2/p10)*drh.drh.st

  d2l.be1.be1  <- -  ( d.n1^2*( -A^2/p11 - A.c^2/p10 - 1/p0)  )  
  d2l.be2.be2  <- -  ( d.n2^2*(B^2*(-1/p11-1/p10))  )
  d2l.be1.be2  <- -  ( d.n1*d.n2*(-A*B/p11+A.c*B/p10)  )
  d2l.be2.rho  <- -  ( -d.n2*d.n1n2*B*(1/p11+1/p10)  )*drh.drh.st 
  d2l.be1.rho  <- -  ( -d.n1*d.n1n2*(A*(1/p11)-A.c*(1/p10))  )*drh.drh.st 
  d2l.rho.rho  <- -  ( -d.n1n2^2*(1/p11+1/p10)   )*drh.drh.st^2 
                                                     
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

  if( ( length(gam1$smooth)==0 && length(gam2$smooth)==0 ) || fp==TRUE){

         list(value=res, gradient=G, hessian=H, l=res, 
              p11=p11, p10=p10, p0=p0, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho) 
  }else{
         S.h <- adiag(matrix(0,gam1$nsdf,gam1$nsdf),
	              S[1:(X1.d2-gam1$nsdf),1:(X1.d2-gam1$nsdf)],
	              matrix(0,gam2$nsdf,gam2$nsdf),
	              S[(X1.d2-(gam1$nsdf-1)):dim(S)[2],(X1.d2-(gam1$nsdf-1)):dim(S)[2]],
                      0)
         S.res <- res
         res <- S.res + 0.5*crossprod(params,S.h)%*%params
         G   <- G + S.h%*%params
         H   <- H + S.h  
         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, 
              p11=p11, p10=p10, p0=p0, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho)  
   }     

}


bprob <- function(params, dat, X1.d2, X2.d2, S=0, gam1, gam2, fp){

  dat1 <- as.matrix(dat[,3:(X1.d2+2)])
  dat2 <- as.matrix(dat[,(X1.d2+3):(X1.d2+X2.d2+2)])
  eta1 <- dat1%*%params[1:X1.d2]
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  corr.st <- params[(X1.d2+X2.d2+1)]
  corr    <- tanh(corr.st)

  y1.y2   <- dat[,1]*dat[,2]
  y1.cy2  <- dat[,1]*(1-dat[,2])
  cy1.y2  <- (1-dat[,1])*dat[,2]
  cy1.cy2 <- (1-dat[,1])*(1-dat[,2])

  d.r <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr^2) )

  A   <- pnorm( (eta2-corr*eta1)*d.r )
  A.c <- pnorm( (eta2-corr*eta1)*d.r , lower.tail = FALSE)
  B   <- pnorm( (eta1-corr*eta2)*d.r )
  B.c <- pnorm( (eta1-corr*eta2)*d.r , lower.tail = FALSE)

  p11 <- pmax( pnorm2( eta1, eta2, corr), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
  p01 <- pmax( pnorm(eta2) - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )

  d.n1   <- dnorm(eta1) 
  d.n2   <- dnorm(eta2) 
  d.n1n2 <- dnorm2(eta1,eta2,corr) 

  l.par <- y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) 

  drh.drh.st   <- 4*exp(2*corr.st)/(exp(2*corr.st)+1)^2
  
  dl.dbe1 <- d.n1*((y1.y2/p11-cy1.y2/p01)*A+(y1.cy2/p10-cy1.cy2/p00)*A.c)
  dl.dbe2 <- d.n2*((y1.y2/p11-y1.cy2/p10)*B+(cy1.y2/p01-cy1.cy2/p00)*B.c)
  dl.drho <- d.n1n2*(y1.y2/p11-y1.cy2/p10-cy1.y2/p01+cy1.cy2/p00)*drh.drh.st

  d2l.be1.be1  <- -(d.n1^2*(A^2*(-1/p11-1/p01)+A.c^2*(-1/p10-1/p00)))      
  d2l.be2.be2  <- -(d.n2^2*(B^2*(-1/p11-1/p10)+B.c^2*(-1/p01-1/p00)))
  d2l.be1.be2  <- -(d.n1*d.n2*(A*B.c/p01-A*B/p11+A.c*B/p10-A.c*B.c/p00))
  d2l.be1.rho  <- -(-d.n1*d.n1n2*(A*(1/p11+1/p01)-A.c*(1/p10+1/p00)))*drh.drh.st 
  d2l.be2.rho  <- -(-d.n2*d.n1n2*(B*(1/p11+1/p10)-B.c*(1/p01+1/p00)))*drh.drh.st 
  d2l.rho.rho  <- -(-d.n1n2^2*(1/p11+1/p01+1/p10+1/p00))*drh.drh.st^2 
                                
  be1.be1 <- t(dat1*c(d2l.be1.be1))%*%dat1
  be2.be2 <- t(dat2*c(d2l.be2.be2))%*%dat2
  be1.be2 <- t(dat1*c(d2l.be1.be2))%*%dat2
  be1.rho <- t(t(rowSums(t(dat1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(dat2*c(d2l.be2.rho)))))

  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
            ) 

  if( ( length(gam1$smooth)==0 && length(gam2$smooth)==0 ) || fp==TRUE){
         res <- -sum(l.par)
         G   <- c( -colSums( c(dl.dbe1)*dat1 ),
                   -colSums( c(dl.dbe2)*dat2 ),
                   -sum( dl.drho )  )
         H   <- H
         list(value=res, gradient=G, hessian=H, l=res, 
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho) 
  }else{
         S.h <- rbind( cbind( cbind(matrix(0,(X1.d2-gam1$nsdf)+gam1$nsdf,gam1$nsdf),rbind(matrix(0,gam1$nsdf,(X1.d2-gam1$nsdf)),S[1:(X1.d2-gam1$nsdf),1:(X1.d2-gam1$nsdf)]))    , matrix(0,dim(be1.be2)[1],dim(be1.be2)[2]) , matrix(0,dim(be1.rho)[1],dim(be1.rho)[2]) ), 
                       cbind( t(matrix(0,dim(be1.be2)[1],dim(be1.be2)[2])) , cbind(matrix(0,length((X1.d2-(gam1$nsdf-1)):dim(S)[2])+gam2$nsdf,gam2$nsdf),rbind(matrix(0,gam2$nsdf,length((X1.d2-(gam1$nsdf-1)):dim(S)[2])),S[(X1.d2-(gam1$nsdf-1)):dim(S)[2],(X1.d2-(gam1$nsdf-1)):dim(S)[2]]))   , matrix(0,dim(be2.rho)[1],dim(be2.rho)[2]) ), 
                       cbind( t(matrix(0,dim(be1.rho)[1],dim(be1.rho)[2])) , t(matrix(0,dim(be2.rho)[1],dim(be2.rho)[2])) , 0 ) ) 
         S.res <- -sum(l.par)
         res <- S.res + (1/2)*(t(params)%*%S.h%*%params)
         G   <- c( -colSums( c(dl.dbe1)*dat1 ),
                   -colSums( c(dl.dbe2)*dat2 ), 
                   -sum(dl.drho) ) + S.h%*%params
         H   <- H + S.h  
         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, 
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho)  
   }     

}


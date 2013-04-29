bprobSS <- function(params, y1, y2, q1, q2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, weights=weights, X1.d2, X2.d2, sp=NULL, qu.mag=NULL, gp1, gp2, fp, l.sp1, l.sp2, K=NULL, n=NULL, N=NULL, cuid=NULL, uidf=NULL, masses=NULL, NGQ=NULL, dat1all=NULL, dat2all=NULL, W=NULL){

  eta1 <- X1%*%params[1:X1.d2]
  eta2 <- X2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  corr.st <- params[(X1.d2+X2.d2+1)]
  corr    <- tanh(corr.st)

  d.r  <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr^2) )

  A   <- pnorm( (eta2-corr*eta1)*d.r )
  A.c <- 1 - A
  B   <- pnorm( (eta1-corr*eta2)*d.r )

  p11 <- pmax( abs(pnorm2( eta1, eta2, cov12=corr)), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
  p0  <- pmax( pnorm(-eta1), 1000*.Machine$double.eps )

  d.n1   <- dnorm(eta1) 
  d.n2   <- dnorm(eta2) 
  d.n1n2 <- dnorm2(eta1,eta2,rho=corr) 

  l.par <- weights*(y1.y2*log(p11)+y1.cy2*log(p10)+cy1*log(p0)) 

  drh.drh.st   <- 4*exp(2*corr.st)/(exp(2*corr.st)+1)^2
  
  dl.dbe1 <- weights*d.n1*( y1.y2/p11*A  + y1.cy2/p10*A.c - cy1/p0 )  
  dl.dbe2 <- weights*d.n2*B*( y1.y2/p11 - y1.cy2/p10)  
  dl.drho <- weights*d.n1n2*(y1.y2/p11 - y1.cy2/p10)*drh.drh.st

  d2l.be1.be1  <- -  weights*( d.n1^2*( -A^2/p11 - A.c^2/p10 - 1/p0)  )  
  d2l.be2.be2  <- -  weights*( d.n2^2*(B^2*(-1/p11-1/p10))  )
  d2l.be1.be2  <- -  weights*( d.n1*d.n2*(-A*B/p11+A.c*B/p10)  )
  d2l.be2.rho  <- -  weights*( -d.n2*d.n1n2*B*(1/p11+1/p10)  )*drh.drh.st 
  d2l.be1.rho  <- -  weights*( -d.n1*d.n1n2*(A*(1/p11)-A.c*(1/p10))  )*drh.drh.st 
  d2l.rho.rho  <- -  weights*( -d.n1n2^2*(1/p11+1/p10)   )*drh.drh.st^2 
                                                     
  be1.be1 <- crossprod(X1*c(d2l.be1.be1),X1)
  be2.be2 <- crossprod(X2*c(d2l.be2.be2),X2)
  be1.be2 <- crossprod(X1*c(d2l.be1.be2),X2)
  be1.rho <- t(t(rowSums(t(X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(X2*c(d2l.be2.rho)))))
  
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
            ) 
            
         res <- -sum(l.par)
         G   <- c( -colSums( c(dl.dbe1)*X1 ),
                   -colSums( c(dl.dbe2)*X2 ),
                   -sum( dl.drho )  )            

if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0

     else{

       S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
       S <- do.call(adiag, lapply(S, unlist)) 

    if(l.sp1!=0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1),
                                           S[1:(X1.d2-gp1),1:(X1.d2-gp1)],
                      			   matrix(0,gp2,gp2),
                      			   S[(X1.d2-(gp1-1)):dim(S)[2],(X1.d2-(gp1-1)):dim(S)[2]],
                      			   0)

    if(l.sp1==0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1), matrix(0,gp2,gp2), S, 0)
    if(l.sp1!=0 && l.sp2==0) S.h <- adiag(matrix(0,gp1,gp1), S, matrix(0,gp2,gp2), 0)


                      
   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
         }

         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
         H   <- H + S.h  

         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, 
              p11=p11, p10=p10, p0=p0, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho)      

}
















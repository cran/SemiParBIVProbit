LM.bpm <- function(formula.eq1,formula.eq2,data=list(),selection=FALSE,FI=FALSE){

G.Eh.E <- function(params,dat,dat1,dat2,X1.d2,X2.d2,gam1,gam2,selection,FI){

  S.h <- 0
  
  if( length(gam1$smooth)!=0 && length(gam2$smooth)!=0 && FI==FALSE){ 
  S <- spS(c(gam1$sp,gam2$sp),gam1,gam2)
  S.h <- adiag(matrix(0,gam1$nsdf,gam1$nsdf),
                      S[1:(X1.d2-gam1$nsdf),1:(X1.d2-gam1$nsdf)],
                      matrix(0,gam2$nsdf,gam2$nsdf),
                      S[(X1.d2-(gam1$nsdf-1)):dim(S)[2],(X1.d2-(gam1$nsdf-1)):dim(S)[2]],
                      0)
  }

  if(FI==TRUE){
  
  eta1 <- dat1%*%params[1:X1.d2]
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  corr <- params[(X1.d2+X2.d2+1)]

  y1.y2   <- dat[,1]*dat[,2]
  y1.cy2  <- dat[,1]*(1-dat[,2])
  
  p11 <- pmax( pnorm2( eta1, eta2, corr), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
  
  d.n1n2 <- dnorm2(eta1,eta2,corr) 
  
  
   if(selection==FALSE){
   
        cy1.y2  <- (1-dat[,1])*dat[,2]
        cy1.cy2 <- (1-dat[,1])*(1-dat[,2])
        p01 <- pmax( pnorm(eta2) - p11, 1000*.Machine$double.eps )
  	p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )      
  
  	dl.drho      <- d.n1n2*(y1.y2/p11-y1.cy2/p10-cy1.y2/p01+cy1.cy2/p00)
  	d2l.rho.rho  <- -(-d.n1n2^2*(1/p11+1/p01+1/p10+1/p00))
  	
   }else{

  	dl.drho      <- d.n1n2*(y1.y2/p11 - y1.cy2/p10) 
        
      	d.r  <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr^2) )
  	A   <- pnorm( (eta2-corr*eta1)*d.r )
  	A.c <- 1 - A
  	B   <- pnorm( (eta1-corr*eta2)*d.r )
  	p0  <- pmax( pnorm(-eta1), 1000*.Machine$double.eps )
  	d.n1   <- dnorm(eta1) 
  	d.n2   <- dnorm(eta2) 

  	d2l.be1.be1  <- -  ( d.n1^2*( -A^2/p11 - A.c^2/p10 - 1/p0)  )  
  	d2l.be2.be2  <- -  ( d.n2^2*(B^2*(-1/p11-1/p10))  )
  	d2l.be1.be2  <- -  ( d.n1*d.n2*(-A*B/p11+A.c*B/p10)  )
  	d2l.be2.rho  <- -  ( -d.n2*d.n1n2*B*(1/p11+1/p10)  )
  	d2l.be1.rho  <- -  ( -d.n1*d.n1n2*(A*(1/p11)-A.c*(1/p10))  )
  	d2l.rho.rho  <- -  ( -d.n1n2^2*(1/p11+1/p10)   )
                          }

  }
  
  if(FI==FALSE){
  
  y1  <- 1
  cy1 <- CY1 <- 0
              
  q1 <- 2*dat[,1]-1
  q2 <- 2*dat[,2]-1
  eta1 <- dat1%*%params[1:X1.d2]
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  corr <- params[(X1.d2+X2.d2+1)]
  corr.sq <- q1*q2*corr
  delta <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr.sq ^2) )

  w1 <- q1*eta1
  w2 <- q2*eta2
  
  g1 <- dnorm(w1)*pnorm( (w2-corr.sq*w1) * delta )
  g2 <- dnorm(w2)*pnorm( (w1-corr.sq*w2) * delta )
  
  PHI2   <- pmax(pnorm2( w1, w2, corr.sq),1000*.Machine$double.eps)
  phi2   <- dnorm2( w1, w2, corr.sq)
  
  v1 <- delta*( w2 - corr.sq*w1 )
  v2 <- delta*( w1 - corr.sq*w2 )


   if(selection==TRUE){
       y1  <- dat[,1]
       cy1 <- (1-y1)  	
       d.n1 <- dnorm(eta1) 
       p0   <- pmax( pnorm(-eta1), 1000*.Machine$double.eps )
       CY1  <- cy1*(-(d.n1/p0)^2 + d.n1*eta1/p0) 
                      }

  dl.drho <- y1*(q1*q2*phi2/PHI2) 
  
  d2l.be1.be1  <- -y1*((-w1*g1 - corr.sq*phi2 - g1^2/PHI2)/PHI2) - CY1
  d2l.be2.be2  <- -y1*((-w2*g2 - corr.sq*phi2 - g2^2/PHI2)/PHI2)
  d2l.be1.be2  <- -y1*((phi2/PHI2 - g1*g2/PHI2^2)*q1*q2)
  d2l.be1.rho  <- -y1*((corr.sq*delta*v1 - w1 - g1/PHI2)*phi2/PHI2*q2)
  d2l.be2.rho  <- -y1*((corr.sq*delta*v2 - w2 - g2/PHI2)*phi2/PHI2*q1)
  d2l.rho.rho  <- -y1*(( phi2/PHI2*( delta^2*corr.sq*( 1 - delta^2*(w1^2+w2^2 - 2*corr.sq*w1*w2) ) + delta^2*w1*w2 - phi2/PHI2 ) ))

  }  
  
  if(FI==FALSE || selection==TRUE){
  
  be1.be1 <- crossprod(dat1*c(d2l.be1.be1),dat1)
  be2.be2 <- crossprod(dat2*c(d2l.be2.be2),dat2)
  be1.be2 <- crossprod(dat1*c(d2l.be1.be2),dat2)
  be1.rho <- t(t(rowSums(t(dat1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(dat2*c(d2l.be2.rho)))))  

  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
             ) + S.h 
                                     
  H.eig <- eigen(H)
  k.e <- sum(as.numeric(H.eig$val<sqrt(.Machine$double.eps)))
     if(k.e!=0){
       ind.e <- (length(H.eig$val)-(k.e-1)):length(H.eig$val)
       min.e <- min(H.eig$val[1:(ind.e[1]-1)])
       for(i in 1:k.e) H.eig$val[ind.e[i]] <- min.e/10^i  
       H.inv <- H.eig$vec%*%tcrossprod(diag(1/H.eig$val),H.eig$vec)      
               }else{
       H.inv <- H.eig$vec%*%tcrossprod(diag(1/H.eig$val),H.eig$vec) 
                    }           
    
  }
  
  if(FI==TRUE && selection==FALSE){
  
    H <- sum(d2l.rho.rho) 
    H.inv <- 1/H
  
  }
                 
  G <- sum( dl.drho )  

L <- list(G=G, V=H, V.inv=H.inv)
L    
}

gam1 <- gam(formula.eq1, binomial(link="probit"), data=data)
 if(selection==FALSE){
	 gam2 <- gam(formula.eq2, binomial(link="probit"), data=data)
	 }else{
	 inde <- gam1$y > 0
	 environment(formula.eq2) <- environment(NULL)
	 gam2  <- gam(formula.eq2, binomial(link="probit"), data=data, subset=inde)
	 environment(gam2$formula) <- environment(gam1$formula)
	 }
	 
params <- c(coef(gam1),coef(gam2),0)

X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]
 if(selection==FALSE){
	X2 <- model.matrix(gam2)
	X2.d2 <- dim(X2)[2]
	dat <- cbind(gam1$y,gam2$y,X1,X2)
	}else{
	X2.d2 <- length(coef(gam2))
	X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
	X2[inde, ] <- model.matrix(gam2)
        y2 <- rep(0,length(inde))
        y2[inde] <- gam2$y
        dat <- cbind(gam1$y,y2,X1,X2)
	}

dat1 <- as.matrix(dat[,3:(X1.d2+2)])
dat2 <- as.matrix(dat[,(X1.d2+3):(X1.d2+X2.d2+2)])

Ev <- G.Eh.E(params,dat,dat1,dat2,X1.d2,X2.d2,gam1,gam2,selection,FI)
if(FI==FALSE || selection==TRUE) ev <- Ev$G^2*Ev$V.inv[(X1.d2+X2.d2+1),(X1.d2+X2.d2+1)]
if(FI==TRUE && selection==FALSE) ev <- Ev$G^2*Ev$V.inv

return(pchisq(ev,1,lower.tail=FALSE))


}


LM.bpm <- function(formula.eq1,formula.eq2,data=list(),selection=FALSE,FI=FALSE){


G.Eh.E <- function(params,y1,y2,X1,X2,X1.d2,X2.d2,gam1,gam2,selection,FI){

  S.h <- 0; l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth) 
  
  if( (l.sp1!=0 || l.sp2!=0) && FI==FALSE){
  
      qu.mag <- S.m(gam1,gam2,l.sp1,l.sp2,K=0,RE=FALSE)

  	if(l.sp1!=0 && l.sp2!=0) sp <- c(gam1$sp,gam2$sp)
  	if(l.sp1==0 && l.sp2!=0) sp <- c(gam2$sp)
  	if(l.sp1!=0 && l.sp2==0) sp <- c(gam1$sp)

      S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
      S <- do.call(adiag, lapply(S, unlist))
      gp1 <- gam1$nsdf; gp2 <- gam2$nsdf 
  
    if(l.sp1!=0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1),
                                           S[1:(X1.d2-gp1),1:(X1.d2-gp1)],
                      			   matrix(0,gp2,gp2),
                      			   S[(X1.d2-(gp1-1)):dim(S)[2],(X1.d2-(gp1-1)):dim(S)[2]],
                      			   0)

    if(l.sp1==0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1), matrix(0,gp2,gp2), S, 0)
    if(l.sp1!=0 && l.sp2==0) S.h <- adiag(matrix(0,gp1,gp1), S, matrix(0,gp2,gp2), 0)

  }


  if(FI==TRUE){
  
  eta1 <- X1%*%params[1:X1.d2]
  eta2 <- X2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  corr <- params[(X1.d2+X2.d2+1)]

  y1.y2   <- y1*y2
  y1.cy2  <- y1*(1-y2)
  
  p11 <- pmax( abs(pbinorm( eta1, eta2, cov12=corr)), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
  
  d.n1n2 <- dbinorm(eta1,eta2,cov12=corr) 
  
  
   if(selection==FALSE){
   
        cy1.y2  <- (1-y1)*y2
        cy1.cy2 <- (1-y1)*(1-y2)
        p01 <- pmax( pnorm(eta2) - p11, 1000*.Machine$double.eps )
  	p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )      
  
  	dl.drho      <- d.n1n2*(y1.y2/p11-y1.cy2/p10-cy1.y2/p01+cy1.cy2/p00)
  	d2l.rho.rho  <- -(-d.n1n2^2*(1/p11+1/p01+1/p10+1/p00))
  	               }
  	               
   if(selection==TRUE){

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

   y1.y2  <- y1*y2
   y1.cy2 <- y1*(1-y2)

  eta1 <- X1%*%params[1:X1.d2]
  eta2 <- X2%*%params[(X1.d2+1):(X1.d2+X2.d2)]
  corr <- params[(X1.d2+X2.d2+1)]
    
  p1 <- pnorm(eta1)
  p2 <- pnorm(eta2)
  
  C.copula <- pbinorm( eta1, eta2, cov12=corr)

  p11 <- pmax( C.copula, 1000*.Machine$double.eps )
  p10 <- pmax( p1 - p11, 1000*.Machine$double.eps )

  d.r    <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-corr^2) )
  d.n1   <- dnorm(eta1) 
  d.n2   <- dnorm(eta2) 
  d.n1n2 <- dbinorm(eta1,eta2, cov12=corr) 

  A   <- pnorm( (eta2-corr*eta1)*d.r )
  B   <- pnorm( (eta1-corr*eta2)*d.r )

  c.copula.be1 <- A  
  c.copula.be2 <- B
  c.copula.theta <- d.n1n2
  c.copula2.be1 <- dnorm((eta2-corr*eta1)*d.r)*-corr*d.r*sqrt(2*pi)/exp(-eta1^2/2)  
  c.copula2.be2 <- dnorm((eta1-corr*eta2)*d.r)*-corr*d.r*sqrt(2*pi)/exp(-eta2^2/2)
  c.copula2.be1be2 <- d.r*exp( -(  corr^2*(eta1^2+eta2^2)-2*corr*eta1*eta2     )/(2*(1-corr^2))          )
  c.copula2.be1th <- -(dnorm((eta2 - corr * eta1)/sqrt(1 - corr^2)) * (eta1/sqrt(1 - 
    corr^2) - (eta2 - corr * eta1) * (0.5 * (2 * corr * (1 - 
    corr^2)^-0.5))/sqrt(1 - corr^2)^2)) 
  c.copula2.be2th <- -(dnorm((eta1 - corr * eta2)/sqrt(1 - corr^2)) * (eta2/sqrt(1 - 
    corr^2) - (eta1 - corr * eta2) * (0.5 * (2 * corr * (1 - 
    corr^2)^-0.5))/sqrt(1 - corr^2)^2))
  bit1.th2 <- 2 * pi * (0.5 * (2 * corr * (1 - corr^2)^-0.5))/(2 * pi * sqrt(1 - 
    corr^2))^2 * exp(-1/(2 * (1 - corr^2)) * (eta1^2 + eta2^2 - 
    2 * corr * eta1 * eta2)) - 1/(2 * pi * sqrt(1 - corr^2)) * 
    (exp(-1/(2 * (1 - corr^2)) * (eta1^2 + eta2^2 - 2 * corr * 
        eta1 * eta2)) * (-1/(2 * (1 - corr^2)) * (2 * eta1 * 
        eta2) + 2 * (2 * corr)/(2 * (1 - corr^2))^2 * (eta1^2 + 
        eta2^2 - 2 * corr * eta1 * eta2))) 



    if(selection==FALSE){

    cy1.y2 <- (1-y1)*y2
    cy1.cy2 <- (1-y1)*(1-y2)


    p01 <- pmax( p2 - p11, 1000*.Machine$double.eps )
    p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )

    dl.drho <- y1.y2*c.copula.theta/p11 + y1.cy2*(-c.copula.theta)/p10 + cy1.y2*(-c.copula.theta)/p01 + cy1.cy2*c.copula.theta/p00  
 
bit1.b1b1 <- c.copula2.be1*(d.n1)^2-c.copula.be1*d.n1*eta1
bit2.b1b1 <- -d.n1*eta1-bit1.b1b1
bit3.b1b1 <- -bit1.b1b1
bit4.b1b1 <- -bit2.b1b1

bit1.b2b2 <- c.copula2.be2*(d.n2)^2-c.copula.be2*d.n2*eta2
bit2.b2b2 <- -bit1.b2b2
bit3.b2b2 <- -d.n2*eta2-bit1.b2b2
bit4.b2b2 <- -bit3.b2b2

bit1.b1b2 <- c.copula2.be1be2 * d.n1 *d.n2
bit2.b1b2 <- -bit1.b1b2
bit3.b1b2 <- -bit1.b1b2
bit4.b1b2 <- bit1.b1b2
       
bit1.b1th <- c.copula2.be1th*d.n1

bit2.b1th <- -bit1.b1th 
bit3.b1th <- -bit1.b1th 
bit4.b1th <- bit1.b1th 


bit1.b2th <- c.copula2.be2th*d.n2
bit2.b2th <- -bit1.b2th 
bit3.b2th <- -bit1.b2th 
bit4.b2th <- bit1.b2th 

bit2.th2 <- -bit1.th2 
bit3.th2 <- -bit1.th2 
bit4.th2 <- bit1.th2 


  d2l.be1.be1  <- -(y1.y2*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2+
                              y1.cy2*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2+
                              cy1.y2*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2+
                              cy1.cy2*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2 )

  d2l.be2.be2  <- -(y1.y2*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              y1.cy2*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2+
                              cy1.y2*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2+
                              cy1.cy2*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2 )

  d2l.be1.be2  <- -(y1.y2*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2+
                              y1.cy2*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2+
                              cy1.y2*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2+
                              cy1.cy2*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2 )

  d2l.be1.rho  <- -(y1.y2*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )

  d2l.be2.rho  <- -(y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )

  d2l.rho.rho  <- -(y1.y2*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                              y1.cy2*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2+
                              cy1.y2*(bit3.th2*p01-(-c.copula.theta)^2)/p01^2+
                              cy1.cy2*(bit4.th2*p00-c.copula.theta^2)/p00^2 )



}


if(selection==TRUE){


  cy1 <- (1-y1)

  p0  <- pmax( pnorm(-eta1), 1000*.Machine$double.eps )

  dl.drho <-  y1.y2*c.copula.theta/p11 + y1.cy2*(-c.copula.theta)/p10 
     
bit1.b1b1 <- c.copula2.be1*(d.n1)^2-c.copula.be1*d.n1*eta1
bit2.b1b1 <- -d.n1*eta1-bit1.b1b1
bit3.b1b1 <- d.n1*(eta1*p0-d.n1)/p0^2
bit1.b2b2 <- c.copula2.be2*(d.n2)^2-c.copula.be2*d.n2*eta2
bit2.b2b2 <- -bit1.b2b2

bit1.b1b2 <- c.copula2.be1be2 * d.n1 *d.n2
bit2.b1b2 <- -bit1.b1b2

bit1.b1th <- c.copula2.be1th*d.n1
bit2.b1th <- -bit1.b1th 


bit1.b2th <- c.copula2.be2th*d.n2
bit2.b2th <- -bit1.b2th 
bit1.th2 <- 2 * pi * (0.5 * (2 * corr * (1 - corr^2)^-0.5))/(2 * pi * sqrt(1 - 
    corr^2))^2 * exp(-1/(2 * (1 - corr^2)) * (eta1^2 + eta2^2 - 
    2 * corr * eta1 * eta2)) - 1/(2 * pi * sqrt(1 - corr^2)) * 
    (exp(-1/(2 * (1 - corr^2)) * (eta1^2 + eta2^2 - 2 * corr * 
        eta1 * eta2)) * (-1/(2 * (1 - corr^2)) * (2 * eta1 * 
        eta2) + 2 * (2 * corr)/(2 * (1 - corr^2))^2 * (eta1^2 + 
        eta2^2 - 2 * corr * eta1 * eta2))) 

bit2.th2 <- -bit1.th2



  d2l.be1.be1  <- -(y1.y2*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2+
                              y1.cy2*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2+
                              cy1*bit3.b1b1 )

  d2l.be2.be2  <- -(y1.y2*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              y1.cy2*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2 )

  d2l.be1.be2  <- -(y1.y2*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2+
                              y1.cy2*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2)

  d2l.be1.rho  <- -(y1.y2*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2 )

  d2l.be2.rho  <- -(y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2 )

  d2l.rho.rho  <- -(y1.y2*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                              y1.cy2*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2 )

}


  }  



  
  if(FI==FALSE || selection==TRUE){
  
  be1.be1 <- crossprod(X1*c(d2l.be1.be1),X1)
  be2.be2 <- crossprod(X2*c(d2l.be2.be2),X2)
  be1.be2 <- crossprod(X1*c(d2l.be1.be2),X2)
  be1.rho <- t(t(rowSums(t(X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(X2*c(d2l.be2.rho)))))  

  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
             ) + S.h 
                                     
  H.eig <- eigen(H,symmetric=TRUE)
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
X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]; y1 <- gam1$y

 if(selection==FALSE){
	 gam2 <- gam(formula.eq2, binomial(link="probit"), data=data)
         X2 <- model.matrix(gam2); X2.d2 <- dim(X2)[2]
         y2 <- gam2$y
	 }else{
	 inde <- y1 > 0
	 gam2 <- eval(substitute(gam(formula.eq2, binomial(link="probit"), data=data, subset=inde),list(inde=inde)))
         X2.d2 <- length(coef(gam2))
	 X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
	 X2[inde, ] <- model.matrix(gam2)
         y2 <- rep(0,length(inde))
         y2[inde] <- gam2$y
	 }
	 
params <- c(coef(gam1),coef(gam2),0)


Ev <- G.Eh.E(params,y1,y2,X1,X2,X1.d2,X2.d2,gam1,gam2,selection,FI)

if(FI==FALSE || selection==TRUE) ev <- Ev$G^2*Ev$V.inv[(X1.d2+X2.d2+1),(X1.d2+X2.d2+1)]
if(FI==TRUE && selection==FALSE) ev <- Ev$G^2*Ev$V.inv

return(pchisq(ev,1,lower.tail=FALSE))

}







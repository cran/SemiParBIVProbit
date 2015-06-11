AT <- function(x, eq, nm.bin, E = TRUE, treat = TRUE, naive = FALSE, ind = NULL, sub.l = 50, delta = FALSE, 
   n.sim = 100, prob.lev = 0.05, hd.plot = FALSE,
   main = "Histogram and Kernel Density of Simulated Average Effects", 
   xlab = "Simulated Average Effects", ...){

etap.noi <- X.int <- X.noi <- eti1 <- eti0 <- etno <- indS <- bs <- ind.excl <- p.int1 <- p.int0 <- d.int1 <- d.int0 <- p.etn <- d.etn <- ass.p <- ass.pst <- C.11 <- C.10 <- sig2 <- peti1s <- peti0s <- sigma2.st <- sigma2s <- eti1s <- eti0s <- d0 <- d1 <- p.etns <- etnos <- etds <- ass.ps <- 1

m2 <- c("N","GU","rGU","LO","LN","WEI","iG","GA")
end <- 0
epsilon <- 0.0000001 # 0.9999999 sqrt(.Machine$double.eps)
max.p   <- 0.9999999
est.ATb <- NA
indD <- list()

if(x$margins[2] != "probit") delta <- FALSE 
if(missing(nm.bin)) stop("You must provide the name of the endogenous variable.")
if(missing(eq) && x$Model=="B" && x$margins[2] != "probit") eq <- 2 
if(x$Model=="B" && x$margins[2] != "probit")                eq <- 2  
if(x$ig[[1]]$response %in% x$ig[[2]]$pred.names ) end <- 1
if(x$ig[[2]]$response %in% x$ig[[1]]$pred.names ) end <- 1
if(x$Model=="BSS" || x$Model=="BPO" || x$Model=="BPO0" || end==0) stop("Calculation of this average treatment effect is valid for recursive models only.")
if(x$gamlssfit == FALSE && naive == TRUE && x$margins[2] %in% m2) stop("You need to fit the naive model to obtain the ATE. Refit the model and set gamlssfit = TRUE.")
if(is.character(nm.bin)==FALSE) stop("nm.bin is not a character!")
if( !is.null(ind) && E == FALSE) stop("ind is not designed to be used when some observations are excluded from the AT's calculation.")

if( !is.null(ind) ){ 

    if(is.logical(ind) == FALSE) stop("ind must be a logical variable.")
    if(length(table(ind))!=2 ) stop("ind must be a logical binary variable.")
    if( length(ind) != x$n ) stop("ind must have the same length as the number of observations used in fitting.")   

}


if(x$margins[2] != "probit" && is.null(ind)){ind <- 1:x$n; ind <- sample(ind, sub.l)}
if(x$margins[2] == "probit" && is.null(ind)) ind <- 1:x$n



if(E == FALSE ) {

 if(x$margins[2] != "probit") ind <- 1:x$n  
 
 if(eq==1) X.int <- as.matrix(x$X1[ind,])
 if(eq==2) X.int <- as.matrix(x$X2[ind,]) 

    if(treat == TRUE)  ind <- as.logical(X.int[, nm.bin]) 
    if(treat == FALSE) ind <- as.logical(X.int[, nm.bin])!=TRUE
                                              
}



########################################################
# Set-up
########################################################


if(naive==FALSE){
	p.rho    <- length(coef(x))
	indD[[1]] <- 1:x$X1.d2 
	indD[[2]] <- x$X1.d2+(1:x$X2.d2)
}


if(eq==1){ X.int <- as.matrix(x$X1[ind,])
    if(naive==FALSE){
           X.noi <- as.matrix(x$X2[ind,])
           ind.int <- indD[[1]]
           ind.noi <- indD[[2]] 
           etap.noi <- x$eta2[ind] 
                    }      
}

if(eq==2){ X.int <- as.matrix(x$X2[ind,])
    if(naive==FALSE){
           X.noi <- as.matrix(x$X1[ind,])
           ind.int <- indD[[2]]
           ind.noi <- indD[[1]] 
           etap.noi  <- x$eta1[ind]
                    } 
}

if(naive==FALSE){              
	   coef.int <- as.numeric(coef(x)[ind.int])
	   coef.noi <- as.numeric(coef(x)[ind.noi])        
}



#                 X.int <- as.matrix(X.int)  # why do I have this?           
#if(naive==FALSE) X.noi <- as.matrix(X.noi)  # same here     
       
                 
d0 <- d1 <- X.int
d0[,nm.bin] <- 0
d1[,nm.bin] <- 1



if(naive==FALSE){
	eti1 <- d1%*%coef.int 
	eti0 <- d0%*%coef.int 
	etno <- etap.noi      
}



#############################################################################
# ATE
#############################################################################

if(naive==FALSE){


##

if(x$margins[2] != "probit"){ 


if( x$margins[2] %in% c("N","GU","rGU","LO") )             { lil <- -Inf;    uil <- Inf}
if( x$margins[2] %in% c("LN","WEI","iG","GA")     )        { lil <- epsilon; uil <- Inf}


###########
# functions
###########

ConExp0 <- function(y2, eti0, sigma2, margin2, p.etn0, ass.p, BivD){   

	pp0 <- distrHsAT(y2, eti0, sigma2, margin2) 
	p2.0   <- pp0$p2
	pdf2.0 <- pp0$pdf2 
	dc0 <- copgHsAT(p1=p.etn0, p2=p2.0, teta=ass.p, BivD=BivD)$c.copula.be2	
	cond0 <- y2*dc0*pdf2.0/p.etn0
	cond0
                                                                   }

ConExp1 <- function(y2, eti1, sigma2, margin2, p.etn0, ass.p, BivD){   

	pp1 <- distrHsAT(y2, eti1, sigma2, margin2)
	p2.1   <- pp1$p2
	pdf2.1 <- pp1$pdf2 	
	dc1 <- copgHsAT(p1=p.etn0, p2=p2.1, teta=ass.p, BivD=BivD)$c.copula.be2	
	cond1 <- y2 * ( (1 - dc1)*pdf2.1/(1-p.etn0))
	cond1
                                                                   }
                                                                   

integr0 <- function(eti0, sigma2, margin2, p.etn0, ass.p, BivD, lil, uil){

  integrate(ConExp0, lower=lil, upper=uil, eti0=eti0, 
            sigma2=sigma2, margin2=margin2, p.etn0=p.etn0, 
            ass.p=ass.p, BivD=BivD)$value
                         
                       }

integr1 <- function(eti1, sigma2, margin2, p.etn0, ass.p, BivD, lil, uil){

  integrate(ConExp1, lower=lil, upper=uil, eti1=eti1,  
            sigma2=sigma2, margin2=margin2, p.etn0=p.etn0, 
            ass.p=ass.p, BivD=BivD)$value
                         
                       }

v.integr0 <- Vectorize(integr0)  
v.integr1 <- Vectorize(integr1) 

###########


ass.p  <- x$theta
sigma2 <- x$sigma2
if( is.null(x$X3) && is.null(x$X4) ) { ass.p <- rep(ass.p, x$n); sigma2 <- rep(sigma2, x$n)     } 

p.etn0  <- pnorm(-etno)
p.etn0 <- pmax(p.etn0, epsilon ) 
p.etn0 <- ifelse(p.etn0 > max.p, max.p, p.etn0)


v0 <- v.integr0(eti0=eti0, sigma2=sigma2[ind], margin2=x$VC$margins[2], p.etn0=p.etn0, 
                ass.p=ass.p[ind], BivD=x$BivD, lil=lil, uil=uil)
v1 <- v.integr1(eti1=eti1, sigma2=sigma2[ind], margin2=x$VC$margins[2], p.etn0=p.etn0, 
                ass.p=ass.p[ind], BivD=x$BivD, lil=lil, uil=uil)

est.ATso <- (v1 - v0)  
   
}  

  

##
if(x$margins[2] == "probit"){ 

p.int1 <- pnorm(eti1)
p.int0 <- pnorm(eti0)

d.int1 <- dnorm(eti1)
d.int0 <- dnorm(eti0)

p.etn  <- pnorm(etno)
d.etn  <- dnorm(etno)


if( is.null(x$X3)  ) { ass.pst <- coef(x)["theta.star"]; ass.p <- x$theta      } 
if( !is.null(x$X3) ) { ass.pst <- x$fit$etad[ind];       ass.p <- x$theta[ind] } 

p.int1 <- pmax(p.int1, epsilon ) 
p.int1 <- ifelse(p.int1 > max.p, max.p, p.int1)
p.int0 <- pmax(p.int0, epsilon )
p.int0 <- ifelse(p.int0 > max.p, max.p, p.int0)

p.etn  <- pmax(p.etn, epsilon )
p.etn <- ifelse(p.etn > max.p, max.p, p.etn)

C.11  <- pmax( BiCDF(p.int1, p.etn, x$nC, ass.p) , epsilon )
C.10  <- p.int0 - pmax( BiCDF(p.int0, p.etn, x$nC, ass.p) , epsilon )
                                 
est.ATso <- ( C.11/p.etn - C.10/(1-p.etn) )  

}   


} # end naive FALSE condition 




##

if(naive == TRUE){

if(x$margins[2] == "probit"){ 

if(eq==1) ngam <- x$gam1
if(eq==2) ngam <- x$gam2

eti1 <- d1%*%coef(ngam) 
eti0 <- d0%*%coef(ngam) 

p.int1 <- pmax(pnorm(eti1), epsilon ) 
p.int1 <- ifelse(p.int1 > max.p,max.p,p.int1) 
p.int0 <- pmax(pnorm(eti0), epsilon ) 
p.int0 <- ifelse(p.int0 > max.p,max.p,p.int0) 

est.ATso <- (p.int1 - p.int0) 

}


if(x$margins[2] != "probit"){ 

parg2 <- x$gamlss$fit$argument[1:x$X2.d2]
eti1  <- d1%*%parg2 
eti0  <- d0%*%parg2
sig2  <- (exp(x$gamlss$fit$sigma2.st) + epsilon)[ind] 

if(x$margins[2] %in% c("N","LO","iG") )   est.ATso <- x$gamlss$fit$argument[nm.bin]  
if(x$margins[2] == "LN")                  est.ATso <- ( exp(eti1 + sig2/2) - exp(eti0 + sig2/2) ) 
if(x$margins[2] == "GU")                  est.ATso <- ( ( eti1 - 0.57722*sqrt(sig2) ) - ( eti0 - 0.57722*sqrt(sig2) ) )   
if(x$margins[2] == "rGU")                 est.ATso <- ( ( eti1 + 0.57722*sqrt(sig2) ) - ( eti0 + 0.57722*sqrt(sig2) ) )                           
if(x$margins[2] == "WEI")                 est.ATso <- (  eti1*gamma((1/sqrt(sig2)) + 1) - eti0*gamma((1/sqrt(sig2)) + 1)  )                          
if(x$margins[2] == "GA")                  est.ATso <- (  exp(eti1) - exp(eti0)   )                           

                 }  



} # end naive condition








###############

est.AT <- mean( est.ATso, na.rm = TRUE )

###############








###############
# DELTA TRUE
###############

if(delta==TRUE){

  if(naive == FALSE){

    if(x$BivD %in% c("N")      )     add.b <- 1/cosh(ass.pst)^2
    if(x$BivD=="F")                  add.b <- 1
    if(x$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) add.b <-  exp(ass.pst)
    if(x$BivD %in% c("C90","C270","J90","J270","G90","G270") ) add.b <- -exp(ass.pst)
   
   dC1 <- copgHs(p1=p.etn,p2=p.int1,teta=ass.p,teta.st=ass.pst,BivD=x$BivD,eta1=etno,eta2=eti1)
   dC0 <- copgHs(p1=p.etn,p2=p.int0,teta=ass.p,teta.st=ass.pst,BivD=x$BivD,eta1=etno,eta2=eti0)

   dATT.noint <- ( (dC1$c.copula.be1*p.etn-C.11)/p.etn^2 + (dC0$c.copula.be1*(1-p.etn)-C.10)/(1-p.etn)^2)*d.etn 
   dATT.noint <- colMeans( (c(dATT.noint)*X.noi) ) 

   dATT.int   <- colMeans( (c( (dC1$c.copula.be2*d.int1)/p.etn )*d1) ) - colMeans( (c( (1-dC0$c.copula.be2)*d.int0/(1-p.etn) )*d0) )

   if( is.null(x$X3) ) dATT.tet   <- mean(((dC1$c.copula.theta/add.b)/p.etn + (dC0$c.copula.theta/add.b)/(1-p.etn)) )  
   
   if( !is.null(x$X3) ) dATT.tet <- colMeans( (c( (dC1$c.copula.theta/add.b)/p.etn + (dC0$c.copula.theta/add.b)/(1-p.etn) )*x$X3[ind,])  )  


   if(eq==2) dATT <- c(dATT.noint,dATT.int,dATT.tet) else dATT <- c(dATT.int,dATT.noint,dATT.tet) 


   var <- bprobgHs(params=coef(x), respvec=x$respvec, VC=x$VC, sp=x$sp, qu.mag=x$qu.mag, AT=TRUE)$hessian


   var.eig <- eigen(var, symmetric=TRUE)                    
   if(min(var.eig$values) < epsilon) var.eig$values[which(var.eig$values < epsilon)] <- 0.0000001
   var <- var.eig$vectors%*%tcrossprod(diag(1/var.eig$values),var.eig$vectors)  

   delta.AT <- sqrt( t(dATT)%*%var%*%dATT )
   
   }
   
   
   if(naive == TRUE){
   
   var <- ngam$Vp   
   dATT <- colMeans(  ( c(dnorm(eti1))*d1 - c(dnorm(eti0))*d0 )  )
   delta.AT <- sqrt( t(dATT)%*%var%*%dATT )
       
   }
   
                      
}





###############
# DELTA FALSE
###############



if(delta==FALSE){     


if(naive == TRUE){ 


if(x$margins[2] == "probit"){ 

  
 bs <- rMVN(n.sim, mean = coef(ngam), sigma=ngam$Vp)
 
 peti1s <- pmax(pnorm(d1%*%t(bs)), epsilon )  
 peti1s <- ifelse(peti1s > max.p,max.p,peti1s)  
 peti0s <- pmax(pnorm(d0%*%t(bs)), epsilon ) 
 peti0s <- ifelse(peti0s > max.p,max.p,peti0s)  
 est.ATb <- colMeans(  (peti1s - peti0s) , na.rm = TRUE    ) 

                             }


if(x$margins[2] != "probit"){ 

 hess <- x$gamlss$fit$hessian                                                                                                     
 He.eig <- eigen(hess, symmetric=TRUE)
 if(min(He.eig$values) < sqrt(.Machine$double.eps)) He.eig$values[which(He.eig$values < sqrt(.Machine$double.eps))] <- 0.0000001
 inv.hess <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values),He.eig$vectors)  

 bs <- rMVN(n.sim, mean = x$gamlss$fit$argument, sigma=inv.hess)
 colnames(bs) <- names(x$gamlss$fit$argument)

 eti1s <- d1%*%t(bs[,1:x$X2.d2])    
 eti0s <- d0%*%t(bs[,1:x$X2.d2])  

 if( !is.null(x$X3) ) sigma2.st <- x$X3[ind,]%*%t(bs[,(x$X2.d2+1):(x$X2.d2+x$X3.d2)]) 
 if(  is.null(x$X3) ) sigma2.st <- bs[,x$X2.d2 + 1]
   
 sigma2.st <- ifelse( sigma2.st > 20, 20, sigma2.st )  
 sigma2.st <- ifelse( sigma2.st < -17, -17, sigma2.st ) 
 sigma2s <- exp(sigma2.st)
               
if(x$margins[2] %in% c("N","LO","iG") )   est.ATb <- bs[,nm.bin] 
if(x$margins[2] == "LN")                  est.ATb <- colMeans(   exp(eti1s + sigma2s/2) - exp(eti0s + sigma2s/2)                         , na.rm = TRUE    )
if(x$margins[2] == "GU")                  est.ATb <- colMeans(    ( eti1s - 0.57722*sqrt(sigma2s) ) - ( eti0s - 0.57722*sqrt(sigma2s) )  , na.rm = TRUE    )   
if(x$margins[2] == "rGU")                 est.ATb <- colMeans(   ( eti1s + 0.57722*sqrt(sigma2s) ) - ( eti0s + 0.57722*sqrt(sigma2s) )   , na.rm = TRUE    )                           
if(x$margins[2] == "WEI")                 est.ATb <- colMeans(    eti1s*gamma(1/sqrt(sigma2s) + 1) - eti0s*gamma(1/sqrt(sigma2s) + 1)    , na.rm = TRUE    )                            
if(x$margins[2] == "GA")                  est.ATb <- colMeans(    exp(eti1s) - exp(eti0s)                                                , na.rm = TRUE    )         
         
                             }

} # end naive loop






if(naive == FALSE){ 

 bs <- rMVN(n.sim, mean = coef(x), sigma=x$Vb) 

 eti1s <- d1%*%t(bs[,ind.int])    
 eti0s <- d0%*%t(bs[,ind.int])   
 etnos <- X.noi%*%t(bs[,ind.noi]) 
 

 if(x$VC$margins[2] %in% m2 ) p.etns  <- pnorm(-etnos) else p.etns  <- pnorm(etnos)
 p.etns  <- pmax(p.etns, epsilon )
 p.etns  <- ifelse(p.etns > max.p,max.p,p.etns) 
 
 
 if(x$VC$margins[2]=="probit"){
   
   if( !is.null(x$X3) ) etds <- x$X3[ind,]%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)])
   if(  is.null(x$X3) ) etds <- bs[,p.rho]
   
   }
   
   if(x$VC$margins[2] %in% m2 ){
   
   if( !is.null(x$X3) ) sigma2.st <- x$X3[ind,]%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]) 
   if(  is.null(x$X3) ) sigma2.st <- bs[,p.rho-1]
   
   if( !is.null(x$X4) ) etds <- x$X4[ind,]%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2 + 1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)])
   if(  is.null(x$X4) ) etds <- bs[,p.rho]  
   
   sigma2.st <- ifelse( sigma2.st > 20, 20, sigma2.st )  
   sigma2.st <- ifelse( sigma2.st < -17, -17, sigma2.st ) 
   sigma2s <- exp(sigma2.st)
  
  } 
 
 
 
 
    if(x$BivD=="F")                   ass.ps <- etds + epsilon
    if(x$BivD %in% c("N"))           {ass.ps <- tanh(etds); ass.ps <- ifelse(ass.ps < -max.p, -max.p, ass.ps) 
                                                            ass.ps <- ifelse(ass.ps >  max.p,  max.p, ass.ps) }
   
    if(x$BivD %in% c("C0", "C180") )  ass.ps <-   exp(etds) + epsilon
    if(x$BivD %in% c("C90","C270") )  ass.ps <- -(exp(etds) + epsilon)
 
    if(x$BivD %in% c("J0", "J180") )  ass.ps <-    1 + exp(etds) + epsilon
    if(x$BivD %in% c("J90","J270") )  ass.ps <- -( 1 + exp(etds) + epsilon)
  
    if(x$BivD %in% c("G0", "G180") )  ass.ps <-    1 + exp(etds)
    if(x$BivD %in% c("G90","G270") )  ass.ps <- -( 1 + exp(etds) )
    
    ass.ps <- ifelse(ass.ps == Inf ,  8.218407e+307, ass.ps )
    ass.ps <- ifelse(ass.ps == -Inf, -8.218407e+307, ass.ps )






if(x$margins[2] != "probit"){                     

if( is.null(x$X3) && is.null(x$X4) ) { ass.ps <- matrix(ass.ps,  ncol = n.sim, nrow=dim(eti0s)[1],byrow = TRUE)
                                      sigma2s <- matrix(sigma2s, ncol = n.sim, nrow=dim(eti0s)[1],byrow = TRUE)  } 

for(i in 1:n.sim){ 

v0 <- v.integr0(eti0=eti0s[,i], sigma2=sigma2s[,i], margin2=x$VC$margins[2], p.etn0=p.etns[,i], ass.p=ass.ps[,i], BivD=x$BivD, lil=lil, uil=uil)
v1 <- v.integr1(eti1=eti1s[,i], sigma2=sigma2s[,i], margin2=x$VC$margins[2], p.etn0=p.etns[,i], ass.p=ass.ps[,i], BivD=x$BivD, lil=lil, uil=uil)

 est.ATb[i] <- mean( (v1 - v0), na.rm = TRUE   ) 

                 }


}







if(x$margins[2] == "probit"){ 

 
p.int1s <- pnorm(eti1s)
p.int0s <- pnorm(eti0s)

p.int1s <- pmax(p.int1s, epsilon )
p.int1s <- ifelse(p.int1s > max.p,max.p,p.int1s)   
p.int0s <- pmax(p.int0s, epsilon )
p.int0s <- ifelse(p.int0s > max.p,max.p,p.int0s) 


if( is.null(x$X3) ) ass.ps <- t(matrix(ass.ps)) 


for(i in 1:n.sim){ 

 C.11 <- BiCDF(p.int1s[,i], p.etns[,i], x$nC, ass.ps[,i]) 
 C.10 <- p.int0s[,i] - BiCDF(p.int0s[,i], p.etns[,i], x$nC, ass.ps[,i]) 

 est.ATb[i] <- mean(   (C.11/p.etns[,i] - C.10/(1-p.etns[,i])), na.rm = TRUE   )
                  }
                  


}   # end probit true






               
}  # end NAIVE FALSE   





  if(hd.plot == TRUE){
  
  if(x$margins[2] %in% m2) mult <- 1 else mult <- 100
  
  hist(est.ATb*mult, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(est.ATb*mult)$y,hist(est.ATb*mult, plot = FALSE)$density)), ...)
  lines(density(est.ATb*mult))

 }






} # end DELTA condition



           
if(delta==TRUE){esd.AT <- delta.AT*qnorm(prob.lev/2,lower.tail = FALSE) 
                   CIs <- c(est.AT - esd.AT, est.AT + esd.AT)
               }else CIs <- as.numeric(quantile(est.ATb,c(prob.lev/2,1-prob.lev/2),na.rm = TRUE))

res <- c(CIs[1],est.AT,CIs[2])
out <- list(res=res, prob.lev=prob.lev, sim.AT=est.ATb, AT.so = est.ATso, mar2=x$margins[2])
 
rm(etap.noi, X.int, X.noi, eti1, eti0, etno, indS, bs, ind.excl, p.int1, p.int0, d.int1, d.int0,
   p.etn, d.etn, ass.p, ass.pst, C.11, C.10, sig2, peti1s, peti0s, sigma2.st, sigma2s, eti1s, eti0s, d0, d1,
   p.etns, etnos, etds, ass.ps)    
 
class(out) <- "AT"

out

}




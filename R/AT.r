AT <- function(x,eq,nm.bin="",E=TRUE,treat=TRUE,delta=TRUE, sig.lev=0.05,s.meth="svd",n.sim=1000){

if(is.character(nm.bin)==FALSE) stop("nm.bin is not a character!")
if(x$sel==TRUE) stop("Calculation of average treatment effects for sample selection models not implemented.")

Eb.u.noi <- Eb.u.int <- 0; est.ATb <- NA; ind <- list()
p.rho <- length(coef(x))
good <- x$fit$good


if(x$RE==FALSE) { ind[[1]] <- 1:x$X1.d2 
                    ind[[2]] <- x$X1.d2+(1:x$X2.d2)
                    }else{
                    ind[[1]] <- (x$K+1):(x$X1.d2+x$K) 
                    ind[[2]] <- (x$X1.d2+2*x$K+1):(x$X1.d2+x$X2.d2+2*x$K) 
                    }


if(eq==1){ X.int <- x$X1
           X.noi <- x$X2
           ind.int <- ind[[1]]
           ind.noi <- ind[[2]] 
           etap.noi <- x$eta2
           
           if(x$RE==TRUE && x$RE.type=="NP"){ 
                             Eb.u.int <- x$Eb.u1
                             Eb.u.noi <- x$Eb.u2
                             }       
}

if(eq==2){ X.int <- x$X2
           X.noi <- x$X1
           ind.int <- ind[[2]]
           ind.noi <- ind[[1]] 
           etap.noi  <- x$eta1
           
           if(x$RE==TRUE && x$RE.type=="NP"){ 
                             Eb.u.int <- x$Eb.u2
                             Eb.u.noi <- x$Eb.u1
                             }  
}


if(x$RE==TRUE && x$RE.type=="NP"){ X.int <- X.int[,-1]
                                   X.noi <- X.noi[,-1]
                                   nw <- x$T.sv$Wp3/matrix(rep(apply(x$T.sv$Wp3,1,sum),x$K),ncol=x$K)}
                  
coef.int <- as.numeric(coef(x)[ind.int])
coef.noi <- as.numeric(coef(x)[ind.noi])        

X.int <- X.int[good,]       
X.noi <- X.noi[good,]   
                             
d0 <- d1 <- X.int
d0[,nm.bin] <- 0
d1[,nm.bin] <- 1

eti1 <- d1%*%coef.int + Eb.u.int 
eti0 <- d0%*%coef.int + Eb.u.int
etno <- etap.noi      + Eb.u.noi


if(E==TRUE) ind.excl <- rep(TRUE,x$n) else{ 

      if(treat==TRUE) ind.excl <- as.logical(X.int[,nm.bin]) 
                 else ind.excl <- as.logical(X.int[,nm.bin])!=TRUE
      
    } 



if(x$gev.eq1==FALSE && x$gev.eq2==FALSE){

p.int1 <- pnorm(eti1)
p.int0 <- pnorm(eti0)

d.int1 <- dnorm(eti1)
d.int0 <- dnorm(eti0)

p.etn  <- pnorm(etno)
d.etn  <- dnorm(etno)

}



if(x$gev.eq1==TRUE && x$gev.eq2==FALSE){

if(eq==1){

p.int1 <- 1-pgev(-eti1,shape=x$shape1)
p.int0 <- 1-pgev(-eti0,shape=x$shape1)

d.int1 <- dgev(-eti1,shape=x$shape1)
d.int0 <- dgev(-eti0,shape=x$shape1)

p.etn  <- pnorm(etno)
d.etn  <- dnorm(etno)}

if(eq==2){

p.int1 <- pnorm(eti1)
p.int0 <- pnorm(eti0)

d.int1 <- dnorm(eti1)
d.int0 <- dnorm(eti0)

p.etn  <- 1-pgev(-etno,shape=x$shape1)
d.etn  <- dgev(-etno,shape=x$shape1) }

}



if(x$gev.eq1==FALSE && x$gev.eq2==TRUE){

if(eq==2){

p.int1 <- 1-pgev(-eti1,shape=x$shape2)
p.int0 <- 1-pgev(-eti0,shape=x$shape2)

d.int1 <- dgev(-eti1,shape=x$shape2)
d.int0 <- dgev(-eti0,shape=x$shape2)

p.etn  <- pnorm(etno)
d.etn  <- dnorm(etno)}

if(eq==1){

p.int1 <- pnorm(eti1)
p.int0 <- pnorm(eti0)

d.int1 <- dnorm(eti1)
d.int0 <- dnorm(eti0)

p.etn  <- 1-pgev(-etno,shape=x$shape2)
d.etn  <- dgev(-etno,shape=x$shape2) }

}



if(x$gev.eq1==TRUE && x$gev.eq2==TRUE){

if(eq==1){

p.int1 <- 1-pgev(-eti1,shape=x$shape1)
p.int0 <- 1-pgev(-eti0,shape=x$shape1)

d.int1 <- dgev(-eti1,shape=x$shape1)
d.int0 <- dgev(-eti0,shape=x$shape1)

p.etn  <- 1-pgev(-etno,shape=x$shape2)
d.etn  <- dgev(-etno,shape=x$shape2)}

if(eq==2){

p.int1 <- 1-pgev(-eti1,shape=x$shape2)
p.int0 <- 1-pgev(-eti0,shape=x$shape2)

d.int1 <- dgev(-eti1,shape=x$shape2)
d.int0 <- dgev(-eti0,shape=x$shape2)

p.etn  <- 1-pgev(-etno,shape=x$shape1)
d.etn  <- dgev(-etno,shape=x$shape1) }

}



if(x$BivD %in% c("N", "T") ) {ass.p <- x$rho; ass.pst <- coef(x)["athrho"] } else{ ass.p <- x$theta; ass.pst <- coef(x)["theta.star"] } 
 
   if(x$BivD=="N") { 
                     pn.int1 <- qnorm(p.int1)
                     pn.int0 <- qnorm(p.int0)
                     pn.etn  <- qnorm(p.etn) 
                     C.11  <- abs(pbinorm(pn.int1,pn.etn,cov12=ass.p))
                     C.10  <- p.int0 - abs(pbinorm(pn.int0,pn.etn,cov12=ass.p))

   }else{ C.11  <- pmax( BiCopCDF(p.int1,p.etn, x$nC, par=ass.p,par2=x$nu) , 1000*.Machine$double.eps )
          C.10  <- p.int0 - pmax( BiCopCDF(p.int0,p.etn, x$nC, par=ass.p, par2=x$nu) , 1000*.Machine$double.eps )
        }

   est.AT <- mean(   (C.11/p.etn - C.10/(1-p.etn))[ind.excl],na.rm = TRUE   )



            
if(delta==TRUE && x$RE==FALSE){ # not sure whether delta method for RE models is possible or it would be quite involved anyway. This is because the estimated linear predictors require quantities which are not calculated directly...

   dC1 <- copgHs(p.etn,p.int1,teta=ass.p,teta.st=ass.pst,x$BivD,x$nC,x$nu)
   dC0 <- copgHs(p.etn,p.int0,teta=ass.p,teta.st=ass.pst,x$BivD,x$nC,x$nu)

   dATT.noint <- ( (dC1$c.copula.be1*p.etn-C.11)/p.etn^2 + (dC0$c.copula.be1*(1-p.etn)-C.10)/(1-p.etn)^2)*d.etn 
   dATT.noint <- colMeans( (c(dATT.noint)*X.noi)[ind.excl,] ) 

   dATT.int   <- colMeans( (c( (dC1$c.copula.be2*d.int1)/p.etn )*d1)[ind.excl,] ) - colMeans( (c( (1-dC0$c.copula.be2)*d.int0/(1-p.etn) )*d0)[ind.excl,])

   dATT.tet   <- mean((dC1$c.copula.theta/p.etn + dC0$c.copula.theta/(1-p.etn))[ind.excl])  

   if(eq==2) dATT <- c(dATT.noint,dATT.int,dATT.tet) else dATT <- c(dATT.int,dATT.noint,dATT.tet) 

   var <- x$Vb   # x$Vb%*%x$HeSh%*%x$Vb

   delta.AT <- sqrt( t(dATT)%*%var%*%dATT )
                       
}else{


if(x$RE==FALSE) {bs <- rmvnorm(n.sim, mean = coef(x), sigma=x$Vb, method=s.meth); Eb.u.ints <- Eb.u.nois <- matrix(0,length(eti1),n.sim)}
              else if(x$RE.type=="NP") bs <- rmvnorm(n.sim, mean = c(coef(x),x$fit$masses[1:(x$K-1)]), sigma=x$Vb, method=s.meth)


if(x$RE==TRUE && x$RE.type=="NP"){

  Eb.u.ints <- Eb.u.nois <- matrix(NA,length(eti1),n.sim)
  
   for(i in 1:n.sim){ 	
                     ind.u1 <- 1:x$K
                     ind.u2 <- (x$K + 1 + x$X1.d2):(2*x$K + x$X1.d2)
	             u1 <- bs[i,ind.u1]
	             u2 <- bs[i,ind.u2]
		     eb.u1 <- apply(t(u1*t(nw)),1,sum)
		     eb.u2 <- apply(t(u2*t(nw)),1,sum) 
		     Eb.u1s <- rep(eb.u1,x$uidf)
		     Eb.u2s <- rep(eb.u2,x$uidf)
		     if(eq==1) {Eb.u.ints[,i] <- Eb.u1s; Eb.u.nois[,i] <- Eb.u2s}  
		     if(eq==2) {Eb.u.ints[,i] <- Eb.u2s; Eb.u.nois[,i] <- Eb.u1s} 		 
                    }
  }


epsilon <- .Machine$double.eps*10^6

for(i in 1:n.sim){ 	
                                        		
 eti1s <- d1%*%bs[i,ind.int]    + Eb.u.ints[,i]
 eti0s <- d0%*%bs[i,ind.int]    + Eb.u.ints[,i]
 etnos <- X.noi%*%bs[i,ind.noi] + Eb.u.nois[,i]


if(x$gev.eq1==FALSE && x$gev.eq2==FALSE){

p.int1s <- pnorm(eti1s)
p.int0s <- pnorm(eti0s)

p.etns  <- pnorm(etnos)

}



if(x$gev.eq1==TRUE && x$gev.eq2==FALSE){

if(eq==1){

p.int1s <- 1-pgev(-eti1s,shape=x$shape1)
p.int0s <- 1-pgev(-eti0s,shape=x$shape1)

p.etns  <- pnorm(etnos)}

if(eq==2){

p.int1s <- pnorm(eti1s)
p.int0s <- pnorm(eti0s)

p.etns  <- 1-pgev(-etnos,shape=x$shape1) }

}



if(x$gev.eq1==FALSE && x$gev.eq2==TRUE){

if(eq==2){

p.int1s <- 1-pgev(-eti1s,shape=x$shape2)
p.int0s <- 1-pgev(-eti0s,shape=x$shape2)

p.etns  <- pnorm(etnos)}

if(eq==1){

p.int1s <- pnorm(eti1s)
p.int0s <- pnorm(eti0s)

p.etns  <- 1-pgev(-etnos,shape=x$shape2) }

}



if(x$gev.eq1==TRUE && x$gev.eq2==TRUE){

if(eq==1){

p.int1s <- 1-pgev(-eti1s,shape=x$shape1)
p.int0s <- 1-pgev(-eti0s,shape=x$shape1)

p.etns  <- 1-pgev(-etnos,shape=x$shape2)}

if(eq==2){

p.int1s <- 1-pgev(-eti1s,shape=x$shape2)
p.int0s <- 1-pgev(-eti0s,shape=x$shape2)

p.etns  <- 1-pgev(-etnos,shape=x$shape1) }

}



 if(x$BivD=="N"){    pn.int1s <- qnorm(p.int1s)
                     pn.int0s <- qnorm(p.int0s)
                     pn.etns  <- qnorm(p.etns) 

est.ATb[i] <- mean(   (abs(pbinorm(pn.int1s,pn.etns,cov12=tanh(bs[i,p.rho])))/p.etns - (p.int0s - abs(pbinorm(pn.int0s,pn.etns,cov12=tanh(bs[i,p.rho]))))/(1-p.etns))[ind.excl],na.rm = TRUE   )

}else{

   if(x$BivD=="F")                   ass.ps <- bs[i,p.rho] + epsilon
   if(x$BivD=="T")                   {ass.ps <- tanh(bs[i,p.rho]); if(ass.ps %in% c(-1,1)) ass.ps <- sign(ass.ps)*0.9999999}
  
   if(x$BivD %in% c("C0", "C180") )  ass.ps <-   exp(bs[i,p.rho]) + epsilon
   if(x$BivD %in% c("C90","C270") )  ass.ps <- -(exp(bs[i,p.rho]) + epsilon)

   if(x$BivD %in% c("J0", "J180") )  ass.ps <-    1 + exp(bs[i,p.rho]) + epsilon
   if(x$BivD %in% c("J90","J270") )  ass.ps <- -( 1 + exp(bs[i,p.rho]) + epsilon)
 
   if(x$BivD %in% c("G0", "G180") )  ass.ps <-    1 + exp(bs[i,p.rho])
   if(x$BivD %in% c("G90","G270") )  ass.ps <- -( 1 + exp(bs[i,p.rho]) )
 
   if(ass.ps=="Inf" || ass.ps=="-Inf" ) ass.ps <- sign(ass.ps)*2.688117e+43

 C.11 <- pmax( BiCopCDF(p.int1s,p.etns, x$nC, par=ass.ps,par2=x$nu) , 1000*.Machine$double.eps )
 C.10 <- p.int0s - pmax( BiCopCDF(p.int0s,p.etns, x$nC, par=ass.ps,par2=x$nu) , 1000*.Machine$double.eps )

 est.ATb[i] <- mean(   (C.11/p.etns - C.10/(1-p.etns))[ind.excl],na.rm = TRUE   )
     }


                 }


}                 

if(delta==TRUE && x$RE==FALSE) {esd.AT <- delta.AT*qnorm(sig.lev/2,lower.tail = FALSE) 
                 CIs <- c(est.AT - esd.AT, est.AT + esd.AT)} else CIs <- as.numeric(quantile(est.ATb,c(sig.lev/2,1-sig.lev/2),na.rm = TRUE))

res <- c(CIs[1],est.AT,CIs[2])
out <- list(res=res, sig.lev=sig.lev)
 
class(out) <- "AT"

out

}




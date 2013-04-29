AT <- function(x,eq,nm.bin="",E=TRUE,treat=TRUE,delta=TRUE, sig.lev=0.05,s.meth="svd",n.sim=1000){

if(is.character(nm.bin)==FALSE) stop("nm.bin is not a character!")
if(x$sel==TRUE) stop("Calculation of average treatment effects for sample selection models not implemented yet.")

Eb.u.noi <- Eb.u.int <- 0; est.ATb <- NA; ind <- list()
p.rho <- length(coef(x))


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

pn.etn <- pnorm(etno)
est.AT <- mean(   (pnorm2(eti1,etno,cov12=x$rho)/pn.etn - pnorm2(eti0,-etno,cov12=-x$rho)/(1-pn.etn))[ind.excl],na.rm = TRUE   )

if(delta==TRUE){

d.r <- 1/sqrt( pmax(10000*.Machine$double.eps, 1-x$rho^2) ); pp <- which(names(coef(x))==nm.bin)
delta.AT <- sqrt( mean( (( (1/pn.etn)*dnorm(eti1)*pnorm( (etno-x$rho*eti1)*d.r ) )^2)[ind.excl]*x$Vb[pp,pp] ) )

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



for(i in 1:n.sim){ 	
                                        		
 eti1s <- d1%*%bs[i,ind.int]    + Eb.u.ints[,i]
 eti0s <- d0%*%bs[i,ind.int]    + Eb.u.ints[,i]
 etnos <- X.noi%*%bs[i,ind.noi] + Eb.u.nois[,i]
 arhos <- bs[i,p.rho] 
 pn.etns <- pnorm(etnos)
 est.ATb[i] <- mean( (pnorm2(eti1s,etnos,cov12=tanh(arhos))/pn.etns - pnorm2(eti0s,-etnos,cov12=-tanh(arhos))/(1-pn.etns))[ind.excl],na.rm = TRUE )
      
                 }
}                 

if(delta==TRUE) {esd.AT <- delta.AT*qnorm(sig.lev/2,lower.tail = FALSE) 
                 CIs <- c(est.AT - esd.AT, est.AT + esd.AT)} else CIs <- quantile(est.ATb,c(sig.lev/2,1-sig.lev/2),na.rm = TRUE)

res <- as.numeric(c(CIs[1],est.AT,CIs[2]))
out <- list(res=res, sig.lev=sig.lev)
 
class(out) <- "AT"

out

}




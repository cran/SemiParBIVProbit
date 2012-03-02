AT <- function(x,eq,nm.bin="",sig.lev=0.05,n.sim=1000,s.meth="svd",E=TRUE,treat=TRUE){

if(is.character(nm.bin)==FALSE) stop("nm.bin is not a character")

Eb.u <- 0; est.ATb <- NA

if(eq==1){ind <- 1:length(x$gam1$coef); coef <- as.numeric(x$fit$argument[ind]); X <- x$X1} 
if(eq==1 && x$npRE==TRUE){ind <- (x$K+1):(x$X1.d2+x$K); coef <- as.numeric(x$fit$argument[ind]); X <- X[,-1]
                          Eb.u <- x$Eb.u1; nw <- x$T.sv$Wp3/matrix(rep(apply(x$T.sv$Wp3,1,sum),3),ncol=x$K)}


if(x$sel==FALSE){

         if(eq==2){ind <- x$X1.d2+(1:length(x$gam2$coef)); coef <- as.numeric(x$fit$argument[ind]); X <- x$X2}
         if(eq==2 && x$npRE==TRUE){ind <- (x$X1.d2+2*x$K+1):(x$X1.d2+x$X2.d2+2*x$K); coef <- as.numeric(x$fit$argument[ind]); 
                                   X <- X[,-1]; Eb.u <- x$Eb.u2; nw <- x$T.sv$Wp3/matrix(rep(apply(x$T.sv$Wp3,1,sum),3),ncol=x$K)}

                }else{

         if(eq==2){ind <- x$X1.d2+(1:length(x$gam2$coef)); coef <- as.numeric(x$fit$argument[ind]); X <- x$X2[x$gam1$y > 0,]}
  }


d0 <- d1 <- X
d0[,nm.bin] <- 0
d1[,nm.bin] <- 1

eta1 <- d1%*%coef + Eb.u 
eta0 <- d0%*%coef + Eb.u 

if(E==TRUE){

        est.AT <- mean( pnorm(eta1) - pnorm(eta0) )

       	   } else {

	 if(treat==TRUE) est.AT <- mean( (pnorm(eta1) - pnorm(eta0))[as.logical(X[,nm.bin])] ) 
                    else est.AT <- mean( (pnorm(eta1) - pnorm(eta0))[as.logical(X[,nm.bin])!=TRUE] )

}


if(x$npRE==FALSE) bs <- rmvnorm(n.sim, mean = x$fit$argument, sigma=x$Vb, method=s.meth)
else bs <- rmvnorm(n.sim, mean = c(x$fit$argument,x$fit$masses[1:(x$K-1)]), sigma=x$Vb, method=s.meth)


if(E==TRUE){

   for(i in 1:n.sim){ 			if(x$npRE==TRUE){
                                        		if(eq==1) ind.u <- 1:x$K else ind.u <- (x$K+x$X1.d2+1):(x$K+x$X1.d2+x$K)
							u <- bs[i,ind.u]
							eb.u <- apply(t(u*t(nw)),1,sum)
							Eb.u <- rep(eb.u,x$uidf)
                                        		}
               est.ATb[i] <- mean( pnorm(d1%*%bs[i,ind] + Eb.u) - pnorm(d0%*%bs[i,ind] + Eb.u) )
                    }

           }else{

if(treat==TRUE) 
   for(i in 1:n.sim){ 			if(x$npRE==TRUE){
                                        		if(eq==1) ind.u <- 1:x$K else ind.u <- (x$K+x$X1.d2+1):(x$K+x$X1.d2+x$K)
							u <- bs[i,ind.u]
							eb.u <- apply(t(u*t(nw)),1,sum)
							Eb.u <- rep(eb.u,x$uidf)
                                        		}
               est.ATb[i] <- mean( (pnorm(d1%*%bs[i,ind] + Eb.u) - pnorm(d0%*%bs[i,ind] + Eb.u))[as.logical(X[,nm.bin])] ) 
                    }


           else for(i in 1:n.sim) {     if(x$npRE==TRUE){
                                        		if(eq==1) ind.u <- 1:x$K else ind.u <- (x$K+x$X1.d2+1):(x$K+x$X1.d2+x$K)
							u <- bs[i,ind.u]
							eb.u <- apply(t(u*t(nw)),1,sum)
							Eb.u <- rep(eb.u,x$uidf)
                                        		}

               est.ATb[i] <- mean( (pnorm(d1%*%bs[i,ind] + Eb.u) - pnorm(d0%*%bs[i,ind] + Eb.u))[as.logical(X[,nm.bin])!=TRUE] )

}

}

CIs <- quantile(est.ATb,c(sig.lev/2,1-sig.lev/2))
res <- as.numeric(cbind(CIs[1],est.AT,CIs[2]))
out <- list(res=res, sig.lev=sig.lev)
 
class(out) <- "AT"

out

}




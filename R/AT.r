AT <- function(x,eq,nm.bin="",sig.lev=0.05,n.sim=1000,s.meth="svd",E=TRUE,treat=TRUE){

if(is.character(nm.bin)==FALSE) stop("nm.bin is not a character")

if(eq==1){ind <- 1:length(x$gam1$coef); coef <- as.numeric(x$fit$argument[ind]); X <- x$X1}

if(x$sel==FALSE){
         if(eq==2){ind <- x$X1.d2+(1:length(x$gam2$coef)); coef <- as.numeric(x$fit$argument[ind]); X <- x$X2}
                }else{
         if(eq==2){ind <- x$X1.d2+(1:length(x$gam2$coef)); coef <- as.numeric(x$fit$argument[ind]); X <- x$X2[x$gam1$y > 0,]}
  }

d0 <- d1 <- X
d0[,nm.bin] <- 0
d1[,nm.bin] <- 1

if(E==TRUE){

	est.AT <- mean( pnorm(d1%*%coef) - pnorm(d0%*%coef) )

       	   } else {

	 if(treat==TRUE) est.AT <- mean( (pnorm(d1%*%coef) - pnorm(d0%*%coef))[as.logical(X[,nm.bin])] ) else est.AT <- mean( (pnorm(d1%*%coef) - pnorm(d0%*%coef))[as.logical(X[,nm.bin])!=TRUE] )

}


bs <- rmvnorm(n.sim, mean = x$fit$argument, sigma=x$Vb, method=s.meth)

est.ATb <- rep(NA,n.sim)


if(E==TRUE){
   for(i in 1:n.sim) est.ATb[i] <- mean(pnorm(d1%*%bs[i,ind])- pnorm(d0%*%bs[i,ind]))
           } else {

if(treat==TRUE) for(i in 1:n.sim) est.ATb[i] <- mean((pnorm(d1%*%bs[i,ind])- pnorm(d0%*%bs[i,ind]))[as.logical(X[,nm.bin])]) else for(i in 1:n.sim) est.ATb[i] <- mean((pnorm(d1%*%bs[i,ind])- pnorm(d0%*%bs[i,ind]))[as.logical(X[,nm.bin])!=TRUE])

}

CIs <- quantile(est.ATb,c(sig.lev/2,1-sig.lev/2))
res <- as.numeric(cbind(CIs[1],est.AT,CIs[2]))
out <- list(res=res, sig.lev=sig.lev)
 
class(out) <- "AT"

out

}










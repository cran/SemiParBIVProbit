AT <- function(x,eq,nm.bin="",sig.lev=0.05,n.sim=1000,s.meth="svd",E=TRUE,treat=TRUE){

if(eq==1){ind <- 1:length(x$gam1$coef); coef <- as.numeric(x$fit$argument[ind]); X <- x$X1}
if(eq==2){ind <- x$X1.d2+(1:length(x$gam2$coef)); coef <- as.numeric(x$fit$argument[ind]); X <- x$X2}

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
res <- round(cbind(CIs[1],est.AT,CIs[2]),3) 

cat("\nAverage Effect with corresponding ",(1-sig.lev)*100,"% 'Confidence' Intervals:\n\n",sep="")
cat(res[2]," (",res[1],",",res[3],")\n\n",sep="")


}










ATE <- function(x,eq,nm.end,sig.lev=0.05,n.sim=1000,s.meth="svd"){

if(eq==1){ind <- 1:length(x$gam1$coef); coef <- as.numeric(x$fit$argument[ind]); d <- data.frame(x$gam1$var.summary)[2,]; l.s <- x$l.sp1}
if(eq==2){ind <- x$X1.d2+(1:length(x$gam2$coef));coef <- as.numeric(x$fit$argument[ind]); d <- data.frame(x$gam2$var.summary)[2,]; l.s <- x$l.sp2}

d <- cbind(inter=1,d)

if(x$l.sp1!=0 && x$l.sp2!=0){
spl.coef <- list()
for(i in 1:l.s) if(eq==1){spl.coef[[i]] <- PredictMat(x$gam1$smooth[[i]], d)}else{spl.coef[[i]] <- PredictMat(x$gam2$smooth[[i]], d)}

if(eq==1) d <- d[,-c((x$gam1$nsdf+1):(x$gam1$nsdf+l.s))]
if(eq==2) d <- d[,-c((x$gam2$nsdf+1):(x$gam2$nsdf+l.s))]

for(i in 1:l.s) d <- cbind(d,spl.coef[[i]])
}

d[nm.end] <- 1
d0 <- d; d0[nm.end] <- 0

d  <- as.numeric(d)
d0 <- as.numeric(d0)

est.ATE <- as.numeric(pnorm(d%*%coef) - pnorm(d0%*%coef))

bs <- rmvnorm(n.sim, mean = x$fit$argument, sigma=x$Vb, method=s.meth)

est.ATEb <- rep(NA,n.sim)

for(i in 1:n.sim) est.ATEb[i] <- pnorm(d%*%bs[i,ind])- pnorm(d0%*%bs[i,ind])

CIs <- quantile(est.ATEb,c(sig.lev/2,1-sig.lev/2))
res <- round(cbind(CIs[1],est.ATE,CIs[2]),3) 

cat("\nAverage Treatment Effect with corresponding ",(1-sig.lev)*100,"% 'Confidence' Intervals:\n\n",sep="")
cat(res[2]," (",res[1],",",res[3],")\n\n",sep="")

}










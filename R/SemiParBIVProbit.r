SemiParBIVProbit <- function(formula.eq1, formula.eq2, data=list(), method="UBRE", pr.tol=1e-10,
                             gamma=1, aut.sp=TRUE, fp=FALSE, start.v=NULL, rinit=1, rmax=100, 
                             fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
        		     control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=.Machine$double.eps^0.5) ){

  #require(VGAM,quietly=TRUE,warn.conflicts=FALSE)
  #require(mgcv,quietly=TRUE,warn.conflicts=FALSE)
  #require(trust,quietly=TRUE,warn.conflicts=FALSE)
  #library(mvtnorm)
  #source("InfCr.r")
  #source("bprob.r")
  #source("spS.r")
  #source("S.m.r")
  #source("plot.SemiParBIVProbit.r")
  #source("summary.SemiParBIVProbit.r")
  #source("print.summary.SemiParBIVProbit.r")
  #source("print.SemiParBIVProbit.r")
  #source("SemiParBIVProbit.r")
  #source("working.comp.r")
  #source("DataGen.r")
  #set.seed(3)
  #dataS <- DataGen(n=500,rho=0.5,coef.bal=FALSE)
  #gamma=1.4
  #formula.eq1 <- y1 ~       s(x1,by=bin.ex) + s(x1,x2)  
  #formula.eq2 <- y2 ~  y1 + bin.ex + s(x1)

  #rinit=1;rmax=100;fterm=sqrt(.Machine$double.eps);mterm=sqrt(.Machine$double.eps)
  #gamma=1;control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=.Machine$double.eps^0.5)
  #data=dataS
  #pr.tol=1e-6
  #aut.sp=TRUE
  #fp=FALSE
  #start.v=NULL

  gam1  <- gam(formula.eq1, binomial(link="probit"), data=data, method="REML")
  gam2  <- gam(formula.eq2, binomial(link="probit"), data=data, method="REML")
  X1    <- model.matrix(gam1); X1.d2 <- dim(X1)[2]
  X2    <- model.matrix(gam2); X2.d2 <- dim(X2)[2]
  l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth); n.e <- 3
  
  dat <- cbind(gam1$y,gam2$y,X1,X2)
  n   <- length(dat[,1])
  
  sp <- c(gam1$sp,gam2$sp)
  if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2),atanh(0.5))

  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2)}
  
  fit  <- trust(bprob, start.v, rinit=rinit, rmax=rmax, dat=dat,
                X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                iterlim=10000, fterm=fterm, mterm=mterm)
                
  conv.sp <- NULL; bs.mgfit <- wor.c <- 0
  
    if(aut.sp==TRUE){

     l.o <- fit$l; l.n <- 0 

      if(l.sp1!=0 && l.sp2!=0){
       j <- 1; conv.sp <- TRUE 
	  while( abs(l.n-l.o)/abs(l.o) > pr.tol ){     

             So <- spS(sp,gam1,gam2)
		 wor.c    <- try(working.comp(fit,X1,X2,X1.d2,X2.d2,n,n.e))
             if(class(wor.c)=="try-error") break
             
             if(method=="UBRE"){
   
                bs.mgfit <- try(magic(y=wor.c$rW.Z,X=wor.c$rW.X,sp=sp,S=qu.mag$Ss,
                                   off=qu.mag$off,rank=qu.mag$rank,n.score=n.e*n,
                                   gcv=FALSE,gamma=gamma,control=control))
                if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                sp <- bs.mgfit$sp
                               }
                               
             if(method=="REML"){   
             
                wY <- wor.c$rW.Z; wX <- wor.c$rW.X
		wS <- vector("list", length(qu.mag$off))
		for(i in 1:length(qu.mag$off)) wS[[i]] <- matrix(0,dim(wX)[2],dim(wX)[2])
		
		bs.d <- NA 
		for(i in 1:length(gam1$smooth)) bs.d[i] <- gam1$smooth[[i]]$df
		j  <- i + 1; jj <- i + length(gam2$smooth)
                for(i in j:jj) bs.d[i] <- gam2$smooth[[i-length(gam1$smooth)]]$df
                bs.d <- bs.d - 1
		  
		for(i in 1:length(qu.mag$off)) 
		        wS[[i]][seq(qu.mag$off[i],qu.mag$off[i]+bs.d[i]),
		                seq(qu.mag$off[i],qu.mag$off[i]+bs.d[i])] <- qu.mag$Ss[[i]] 
		  
		vSn <- ""; for (i in 1:length(qu.mag$off)) vSn[i] <- paste("wS[[",i,"]]",sep="")
		  
		fw <- paste("list(wX=list(",vSn[1],sep="")
		for (i in 2:length(qu.mag$off)) fw <- paste(fw,",",vSn[i],sep="")
		fw <- paste(fw,"))",sep="")
		pP <- eval(parse( text=fw ))
		bs.mgfit <- gam(wY ~ wX - 1, paraPen=pP, method="REML")
		if(class(bs.mgfit)[1]=="try-error") {conv.sp <- FALSE; break} 
                sp <- bs.mgfit$sp
                               }
		             
             S <- spS(sp,gam1,gam2)
             
             l.o <- fit$l
             fit <- try(trust(bprob, fit$argument, rinit=rinit, rmax=rmax, dat=dat, X1.d2=X1.d2, X2.d2=X2.d2, S=S, 
                          gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, iterlim=10000, fterm=fterm, mterm=mterm),silent=TRUE)

              if(class(fit)=="try-error"){ 
               fit  <- trust(bprob, c(coef(gam1),coef(gam2),atanh(0.5)), rinit=rinit, rmax=rmax, dat=dat,
                             X1.d2=X1.d2, X2.d2=X2.d2, S=So, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                             iterlim=10000, fterm=fterm, mterm=mterm)
               conv.sp <- FALSE
               break
              }

             l.n <- fit$l

             if(j>50){
              conv.sp <- FALSE
              break
             }
             j <- j + 1       
           
        }
      }
    }  

  He <- fit$hessian
  He.eig <- eigen(He)
  k.e <- sum(as.numeric(He.eig$val<.Machine$double.eps^0.4))
  
   if(k.e!=0){
    ind.e <- (length(He.eig$val)-(k.e-1)):length(He.eig$val)
    min.e <- min(He.eig$val[1:(ind.e[1]-1)])
    for(i in 1:k.e) He.eig$val[ind.e[i]] <- min.e/10^i  
    Vb <- He.eig$vec%*%diag(1/He.eig$val)%*%t(He.eig$vec)      
   }else{
    Vb <- He.eig$vec%*%diag(1/He.eig$val)%*%t(He.eig$vec) 
    }
             
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){HeSh <- He - fit$S.h; F <- Vb%*%HeSh} else{HeSh <- He; F <- Vb%*%HeSh}      
  t.edf <- sum(diag(F))

L <- list(fit=fit, gam1=gam1, gam2=gam2, sp=sp, X1=X1, X2=X2,
          rho=tanh(fit$argument[length(fit$argument)]), n=n, X1.d2=X1.d2, X2.d2=X2.d2, 
          l.sp1=l.sp1, l.sp2=l.sp2, n.e=n.e, He=He, Vb=Vb, HeSh=HeSh, F=F, 
          t.edf=t.edf, bs.mgfit=bs.mgfit, conv.sp=conv.sp, wor.c=wor.c,
          p11=fit$p11, p10=fit$p10, p01=fit$p01, p00=fit$p00, eta1=fit$eta1, eta2=fit$eta2, dat=dat)

class(L) <- "SemiParBIVProbit"

L

}

SemiParBIVProbit <- function(formula.eq1, formula.eq2, data=list(), gcv=FALSE, selection=FALSE, 
                             iterlimFS=1, iterlimSP=25, pr.tol=1e-6,
                             gamma=1, aut.sp=TRUE, fp=FALSE, start.v=NULL, rinit=1, rmax=100, 
                             fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
        		     control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=.Machine$double.eps^0.5) ){

  gam1 <- gam(formula.eq1, binomial(link="probit"), gamma=gamma, data=data)
  conv.sp <- gam2.1 <- NULL; bs.mgfit <- wor.c <- j.it <- 0

  if(selection==FALSE){

  gam2  <- gam(formula.eq2, binomial(link="probit"), gamma=gamma, data=data)
  X1    <- model.matrix(gam1); X1.d2 <- dim(X1)[2]
  X2    <- model.matrix(gam2); X2.d2 <- dim(X2)[2]
  l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth)
  dat <- cbind(gam1$y,gam2$y,X1,X2); n <- length(dat[,1])

  i.rho <- 0.5; names(i.rho) <- "rho" 
  if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2),atanh(i.rho))

  func.opt <- bprob
  sp <- c(gam1$sp,gam2$sp)
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2)}

                     }else{

  inde <- gam1$y > 0
  environment(formula.eq2) <- environment(NULL)
  gam2  <- gam(formula.eq2, binomial(link="probit"), gamma=gamma, data=data, subset=inde)
  environment(gam2$formula) <- environment(gam1$formula)
  X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]
  X2.d2 <- length(coef(gam2)); X2 <- matrix(0,length(inde),X2.d2)
  X2[inde,] <- model.matrix(gam2)
  y2 <- rep(0,length(inde)); y2[inde] <- gam2$y
  l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth)
  dat <- cbind(gam1$y,y2,X1,X2); n <- length(dat[,1])

         if(is.null(start.v)){ p.g1 <- predict(gam1)
  	 		       imr <- dnorm(p.g1)/pnorm(p.g1)
  		               formula.eq2.1 <- update.formula(formula.eq2, ~. + imr)
  			       environment(formula.eq2.1) <- environment(NULL)
                               ols <- gam(formula.eq2.1, data=data, gamma=gamma, subset=inde) 
                               rhi <- coef(ols)["imr"]/sqrt(ols$sig2)
                               ta2 <- (1 + rhi^2*imr*(-p.g1-imr))[inde]

                               M <- ols$model
                               vSn <- names(M); fw <-  paste(vSn[1],"~",vSn[2],sep="") 
                               for (i in 3:length(vSn)) fw <- paste(fw,"+",vSn[i],sep="")
                               fw <- as.formula(fw)
                               MM <- as.data.frame(model.matrix(fw,M)[,-1])
                               data.r <- MM/sqrt( pmax(10000*.Machine$double.eps, ta2) ) 
                               data.r <- as.data.frame(cbind(ols$model[,1],data.r)); names(data.r)[1] <- names(M)[1]

                               l <- as.matrix(data.r[,2:(ols$smooth[[1]]$first.para-1)])
                               fw <-  paste(names(data.r)[1],"~ l",sep="")
                               for(i in 1:length(ols$smooth)) fw <- paste(fw,"+",ols$smooth[[i]]$label,sep="")
                               fw <- as.formula(fw)

	                       gam2.1 <- gam(fw, binomial(link="probit"), gamma=gamma, data=data.r)
                               names(gam2.1$coefficients)[2:(ols$smooth[[1]]$first.para-1)] <- names(ols$coefficients)[2:(ols$smooth[[1]]$first.para-1)]

                               rho.c <- coef(gam2.1)["imr"]  
                               rho <- ifelse( abs(rho.c) > 0.99, sign(rho.c)*0.95, rho.c ); names(rho) <- "rho"
  			       start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], atanh(rho)  )

  			     }
  func.opt <- bprobSS
  sp <- c(gam1$sp,gam2.1$sp)
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2)}

                          }


    fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, dat=dat,
                  X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                  fterm=fterm, mterm=mterm)   

                                                     
    if(aut.sp==TRUE){

     l.o <- fit$l; l.n <- 0 

      if(l.sp1!=0 && l.sp2!=0){
      
       j.it <- 1; conv.sp <- TRUE 

	  while( abs(l.n-l.o)/abs(l.o) > pr.tol ){     

             So <- spS(sp,gam1,gam2); coefo <- fit$argument   
  
		 wor.c <- try(working.comp(fit,X1,X2,X1.d2,X2.d2,n))
                 if(class(wor.c)=="try-error") break
             
                	bs.mgfit <- try(magic(y=wor.c$rW.Z,X=wor.c$rW.X,sp=sp,S=qu.mag$Ss,
                        	           off=qu.mag$off,rank=qu.mag$rank,n.score=3*n,
                                	   gcv=gcv,gamma=gamma,control=control))
                	if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                	sp <- bs.mgfit$sp
                               
             S <- spS(sp,gam1,gam2)
             
             l.o <- fit$l
             
             fit <- try(trust(func.opt, fit$argument, rinit=rinit, rmax=rmax, dat=dat, X1.d2=X1.d2, X2.d2=X2.d2, S=S, 
                              gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, iterlim=iterlimFS, fterm=fterm, mterm=mterm),silent=TRUE)

              if(class(fit)=="try-error"){ 
               fit  <- trust(func.opt, coefo, rinit=rinit, rmax=rmax, dat=dat,
                             X1.d2=X1.d2, X2.d2=X2.d2, S=So, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                             iterlim=iterlimFS, fterm=fterm, mterm=mterm)
               conv.sp <- FALSE; break
              }
              
             l.n <- fit$l

             if(j.it>iterlimSP){
              conv.sp <- FALSE; break
             }
             j.it <- j.it + 1      
           
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

L <- list(fit=fit, gam1=gam1, gam2=gam2, gam2.1=gam2.1, sp=sp, iter.sp=j.it,
          rho=tanh(fit$argument[length(fit$argument)]), n=n, X1=X1, X2=X2, X1.d2=X1.d2, X2.d2=X2.d2, 
          l.sp1=l.sp1, l.sp2=l.sp2, He=He, HeSh=HeSh, Vb=Vb, F=F, 
          t.edf=t.edf, bs.mgfit=bs.mgfit, conv.sp=conv.sp, wor.c=wor.c,
          p11=fit$p11, p10=fit$p10, p01=fit$p01, p00=fit$p00, p0=fit$p0, eta1=fit$eta1, eta2=fit$eta2, 
          dat=dat,sel=selection)

class(L) <- "SemiParBIVProbit"

L

}

SemiParBIVProbit <- function(formula.eq1, formula.eq2, data=list(), selection=FALSE, 
                             iterlimSP=50, pr.tol=1e-6,
                             gamma=1, aut.sp=TRUE, fp=FALSE, start.v=NULL, rinit=1, rmax=100, 
                             fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
        		     control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=sqrt(.Machine$double.eps)),
                             npRE=FALSE, K=3, id=NULL){
                             
  conv.sp <- gam2.1 <- bs.mgfit <- wor.c <- j.it <- count.npRE <- N <- cuid <- uidf <- masses <- logL.RE <- NULL

  gam1 <- gam(formula.eq1, binomial(link="probit"), gamma=gamma, data=data)
  
  if(selection==FALSE){

  gam2  <- gam(formula.eq2, binomial(link="probit"), gamma=gamma, data=data)
  X1    <- model.matrix(gam1); X1.d2 <- dim(X1)[2]
  X2    <- model.matrix(gam2); X2.d2 <- dim(X2)[2]
  l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth)
  dat <- cbind(gam1$y,gam2$y,X1,X2); n <- length(dat[,1])

  i.rho <- 0.5; names(i.rho) <- "rho" 
  if(is.null(start.v) && npRE==FALSE) start.v <- c(coef(gam1),coef(gam2),atanh(i.rho))

  func.opt <- bprob
  sp <- c(gam1$sp,gam2$sp)


      if(npRE==TRUE){ 
      
           if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2)}
      
           if(is.null(start.v)){ 
           
           start.v <- c(coef(gam1),coef(gam2),atanh(i.rho))
           fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, dat=dat,
                         X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                         fterm=fterm, mterm=mterm, iterlim = 1e+4,
                         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses) 
           start.v <- fit$argument
                                }
                               

      	dat <- dat[,-c(3,(3+X1.d2))]
        X1.d2 <- X1.d2-1; X2.d2 <- X2.d2-1
        
        	start.v <- start.v[-c(1,X1.d2+2)]
         	start.v <- c(sqrt(2)*gauss.quad(K,"hermite")$nodes,start.v[1:X1.d2],
                      	     sqrt(2)*gauss.quad(K,"hermite")$nodes,start.v[(X1.d2+1):(X2.d2+X1.d2)],
                      	     sign(start.v[X2.d2+X1.d2+1])*i.rho)
                for(i in 1:K) {names(start.v)[i] <- paste("re",i,sep="")
                               names(start.v)[(X1.d2+K+i)] <- paste("re",i,sep="")
                              }	
                nsv <- names(start.v)              
                              
        uid <- unique(id) 
      	N <- length(uid) # how many individuals in the panel
      	uidf <- array(0,N) 
      	for (j in 1:N) uidf[j] <- sum(id==unique(id)[j]) # number of observations in each cluster
      	cuid <- c(0,cumsum(uidf)) # cumulative function used for internal calculations

      	masses <- rep(1/K,K) 
      	func.opt <- bprobNP

        T.sv <- func.opt(start.v, dat=dat, 
                         X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, 
                         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)
                         
          count.npRE <- 0 

      	  while( max(abs(T.sv$gradient)) > pr.tol*10000){
	    start.v <- start.v - ginv(T.sv$hessian)%*%T.sv$gradient
	    masses <- T.sv$masses
	    T.sv <- func.opt(start.v, dat=dat, 
                             X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, 
                             K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)
	    count.npRE <- count.npRE + 1
	    #print(T.sv$W)
                                  }
            start.v <- as.vector(start.v); names(start.v) <- nsv                       
                                  
                    } 

                     }else{

  inde <- gam1$y > 0
  environment(formula.eq2) <- environment(NULL)
  gam2  <- gam(formula.eq2, binomial(link="probit"), gamma=gamma, data=data, subset=inde)
  environment(gam2$formula) <- environment(gam1$formula)
  X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]; X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2)
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
                               if(dim(M)[2]>2) for (i in 3:length(vSn)) fw <- paste(fw,"+",vSn[i],sep="")
                               fw <- as.formula(fw)
                               MM <- as.data.frame(model.matrix(fw,M)[,-1])
                               if(dim(M)[2]==2) names(MM) <- vSn[2]
                               data.r <- MM/sqrt( pmax(10000*.Machine$double.eps, ta2) ) 
                               data.r <- as.data.frame(cbind(ols$model[,1],data.r)); names(data.r)[1] <- names(M)[1]
                        
                               if(l.sp1!=0 && l.sp2!=0 && fp==FALSE) l <- as.matrix(data.r[,2:(ols$smooth[[1]]$first.para-1)]) else l <- as.matrix(data.r[,-1])
                               fw <-  paste(names(data.r)[1],"~ l",sep="")
                               if(l.sp1!=0 && l.sp2!=0 && fp==FALSE) for(i in 1:length(ols$smooth)) fw <- paste(fw,"+",ols$smooth[[i]]$label,sep="")
                               fw <- as.formula(fw)
	                       gam2.1 <- suppressWarnings(try(gam(fw, binomial(link="probit"), gamma=gamma, data=data.r),silent=TRUE))

                               if(class(gam2.1)[1]=="try-error"){
                                                                 names(rhi) <- "rho" 
                                                                 start.v <- c(coef(gam1),coef(gam2),atanh(rhi))
                                                                 gam2.1 <- gam2 
                                                                } else {

                               environment(gam2.1$formula) <- environment(gam1$formula)
                               if(l.sp1!=0 && l.sp2!=0 && fp==FALSE) names(gam2.1$coefficients)[1:(ols$smooth[[1]]$first.para-1)] <- names(ols$coefficients)[1:(ols$smooth[[1]]$first.para-1)] else names(gam2.1$coefficients) <- names(ols$coefficients) 
                               rho.c <- coef(gam2.1)["imr"]  
                               rho <- ifelse( abs(rho.c) > 0.99, sign(rho.c)*0.95, rho.c ); names(rho) <- "rho"
  			       start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], atanh(rho)  )
                                                                       }

  			     }
  func.opt <- bprobSS
  sp <- c(gam1$sp,gam2.1$sp)
  

                          }

    if(l.sp1!=0 && l.sp2!=0 && fp==FALSE && npRE!=TRUE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2)}
    
    
    fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, dat=dat,
                  X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                  fterm=fterm, mterm=mterm, iterlim = 1e+4, 
                  K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)  

    iter.if <- fit$iterations    
                                             
    if(aut.sp==TRUE){

      if(l.sp1!=0 && l.sp2!=0){
      
       j.it <- stoprule.SP <- 1; conv.sp <- TRUE 

	  while( stoprule.SP > pr.tol ){ 

             So <- spS(sp,gam1,gam2); coefo <- fit$argument   
  
		 wor.c <- try(working.comp(fit,X1,X2,X1.d2,X2.d2,n))
                 if(class(wor.c)=="try-error") break
             
                	bs.mgfit <- try(magic(y=wor.c$rW.Z,X=wor.c$rW.X,sp=sp,S=qu.mag$Ss,
                        	           off=qu.mag$off,rank=qu.mag$rank,n.score=3*n,
                                	   gcv=FALSE,gamma=gamma,control=control))
                	if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                	sp <- bs.mgfit$sp
                               
             S <- spS(sp,gam1,gam2)
             o.ests <- c(fit$argument)
            
             fit <- try(trust(func.opt, fit$argument, rinit=rinit, rmax=rmax, dat=dat, X1.d2=X1.d2, X2.d2=X2.d2, S=S, 
                              gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, iterlim=1e+4, fterm=fterm, mterm=mterm),silent=TRUE)

              if(class(fit)=="try-error"){ 
               fit  <- trust(func.opt, coefo, rinit=rinit, rmax=rmax, dat=dat,
                             X1.d2=X1.d2, X2.d2=X2.d2, S=So, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                             iterlim=1e+4, fterm=fterm, mterm=mterm)
               conv.sp <- FALSE; break
                                          } 
              
             n.ests <- c(fit$argument)

             if(j.it>iterlimSP){
              conv.sp <- FALSE; break
                                }
             stoprule.SP <- max(abs(o.ests-n.ests))
             j.it <- j.it + 1      
           
          }
      }
    }
    
    if(selection==FALSE && npRE==TRUE){
    
        func.opt <- bprobNP.H
        T.sv <- func.opt(fit$argument, dat=dat,
                         X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp,
                         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)

        dat1 <- function(u) cbind(matrix(replace(rep(0,K),u,1),ncol=K,nrow=n,byrow=TRUE),T.sv$dat1p)
        dat2 <- function(u) cbind(matrix(replace(rep(0,K),u,1),ncol=K,nrow=n,byrow=TRUE),T.sv$dat2p)

        npar <- K + X1.d2 + K + X2.d2 + 1 + K - 1
        H.cor <- G.cor <- matrix(0,nrow=npar,ncol=npar) #check if we can do it in a faster way

        for (i in 1:N){
           h <- matrix(0,nrow=npar,ncol=1)
           for (l in 1:K){

              c1 <- T.sv$dl.dbe1[(cuid[i]+1):(cuid[i+1]),l]*dat1(l)[(cuid[i]+1):(cuid[i+1]),]
              if (uidf[i]>1) c1<-apply(c1,2,sum)

              c2 <- T.sv$dl.dbe2[(cuid[i]+1):(cuid[i+1]),l]*dat2(l)[(cuid[i]+1):(cuid[i+1]),]
              if (uidf[i]>1) c2<-apply(c2,2,sum)

              c3 <- sum(T.sv$dl.drho[(cuid[i]+1):(cuid[i+1]),l])

              if (l < K) {c4 <- array(0,K-1); c4[l] <- (1/T.sv$masses[l])}else{ c4 <- array((-1/T.sv$masses[K]),K-1)}

              g <- c(c1,c2,c3,c4)
              g <- as.matrix(g)

              H.cor <- H.cor + T.sv$W[i,l]*(g%*%t(g))

              h <- h + T.sv$W[i,l]*g
           }
           G.cor <- G.cor + h%*%t(h)
        }

        I.prob1 <- apply(t((1/T.sv$masses^2)*t(T.sv$W)),2,sum)
        I.prob <- diag(I.prob1[1:K-1])+matrix(I.prob1[K],K-1,K-1)

        He <- adiag(T.sv$hessian,I.prob) - H.cor + G.cor

        #log-Likelihood
        Wp1 <- exp(T.sv$l.par)
        Wp2 <- matrix(0,nrow=N,ncol=K)
        for (i in 1:N) if (uidf[i]>1) {Wp2[i,] <- apply(Wp1[(cuid[i]+1):(cuid[i+1]),],2,prod)
                                 }else{Wp2[i,] <- Wp1[(cuid[i]+1):(cuid[i+1]),]}
        Wp3 <- t(masses*t(Wp2))
        logL.RE <- sum(log(apply(Wp3,1,sum)))
        
    }else{
    
    He <- fit$hessian}
    
    He.eig <- eigen(He)
    k.e <- sum(as.numeric(He.eig$val<sqrt(.Machine$double.eps)))
    
     if(k.e!=0){
      ind.e <- (length(He.eig$val)-(k.e-1)):length(He.eig$val)
      min.e <- min(He.eig$val[1:(ind.e[1]-1)])
      for(i in 1:k.e) He.eig$val[ind.e[i]] <- min.e/10^i  
      Vb <- He.eig$vec%*%diag(1/He.eig$val)%*%t(He.eig$vec)      
     }else{
      Vb <- He.eig$vec%*%diag(1/He.eig$val)%*%t(He.eig$vec) 
    }
                                              
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){ if(npRE==TRUE) fit$S.h <- adiag(fit$S.h,matrix(0,K-1,K-1)) 
                                             HeSh <- He - fit$S.h; F <- Vb%*%HeSh
                                       }else{HeSh <- He; F <- Vb%*%HeSh}      
  t.edf <- sum(diag(F))

L <- list(fit=fit, gam1=gam1, gam2=gam2, gam2.1=gam2.1, sp=sp, iter.sp=j.it, iter.if=iter.if,
          rho=tanh(fit$argument[length(fit$argument)]), n=n, n.sel=length(gam2$y), 
          X1=X1, X2=X2, X1.d2=X1.d2, X2.d2=X2.d2, 
          l.sp1=l.sp1, l.sp2=l.sp2, He=He, HeSh=HeSh, Vb=Vb, F=F, 
          t.edf=t.edf, bs.mgfit=bs.mgfit, conv.sp=conv.sp, wor.c=wor.c,
          p11=fit$p11, p10=fit$p10, p01=fit$p01, p00=fit$p00, p0=fit$p0, eta1=fit$eta1, eta2=fit$eta2, 
          dat=dat,sel=selection,masses=fit$masses,K=K,iter.npRE=count.npRE, npRE=npRE, logL.RE=logL.RE)

class(L) <- "SemiParBIVProbit"

L

}

SemiParBIVProbit <- function(formula.eq1, formula.eq2, data=list(), p.weights=NULL, selection=FALSE, H=FALSE, 
                             iterlimSP=50, pr.tol=1e-6,
                             gamma=1, aut.sp=TRUE, fp=FALSE, start.v=NULL, rinit=1, rmax=100, 
                             fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
        		     control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=sqrt(.Machine$double.eps)),
                             npRE=FALSE, K=3, id=NULL, e.npRE=TRUE, e.npREsp=TRUE){
                             
  conv.sp <- startvSS <- bs.mgfit <- wor.c <- j.it <- count.npRE <- count.npREsp <- N <- cuid <- uidf <- masses <- logL.RE <- eb.u1 <- eb.u2 <- Eb.u1 <- Eb.u2 <- uidf <- T.sv <- NULL
             
  if(!is.null(p.weights)){
            p.weights <- as.vector(p.weights)
            if (!is.numeric(p.weights)) stop("prior weights must be a numeric vector")
            if (any(p.weights < 0))     stop("negative prior weights not allowed")
      }
        
  environment(formula.eq1) <- environment(NULL)                                                      
  gam1 <- gam(formula.eq1, binomial(link="probit"), gamma=gamma, weights=p.weights, data=data)
  environment(gam1$formula) <- environment(formula.eq2)  

  if(is.null(p.weights)) p.weights <- rep(1,length(gam1$y))
  
  if(selection==FALSE){

  environment(formula.eq2) <- environment(NULL) 
  gam2  <- gam(formula.eq2, binomial(link="probit"), gamma=gamma, weights=p.weights, data=data)
  environment(gam2$formula) <- environment(gam1$formula) 

  X1    <- model.matrix(gam1); X1.d2 <- dim(X1)[2]
  X2    <- model.matrix(gam2); X2.d2 <- dim(X2)[2]
  l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth)
  dat <- cbind(gam1$y,gam2$y,X1,X2); n <- length(dat[,1])

  i.rho <- 0.5; names(i.rho) <- "rho" 
  if(is.null(start.v) && npRE==FALSE) start.v <- c(coef(gam1),coef(gam2),atanh(i.rho))

  if(H==TRUE && l.sp1==0 && l.sp2==0) func.opt <- bprobH else func.opt <- bprob
  
  sp <- c(gam1$sp,gam2$sp)

  dat1 <- as.matrix(dat[,3:(X1.d2+2)])
  dat2 <- as.matrix(dat[,(X1.d2+3):(X1.d2+X2.d2+2)])

      if(npRE==TRUE){ 

           if(dim(dat)[1]!=length(id)) stop("The length of your id does not match the number of rows of your data frame")
      
           if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2,K,npRE=FALSE)}
      
           if(is.null(start.v)){ 
           
           start.v <- c(coef(gam1),coef(gam2),atanh(i.rho))
           fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, dat=dat, dat1=dat1, dat2=dat2,
                         X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                         fterm=fterm, mterm=mterm, iterlim = 1e+4, p.weights=p.weights,
                         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses) 
           start.v <- fit$argument
                                }
                               
      	dat <- dat[,-c(3,(3+X1.d2))]
        X1.d2 <- X1.d2-1; X2.d2 <- X2.d2-1
        
        	start.v <- start.v[-c(1,X1.d2+2)]
         	start.v <- c(sqrt(2)*gauss.quad(K,"hermite")$nodes,start.v[1:X1.d2],
                      	     sqrt(2)*gauss.quad(K,"hermite")$nodes,start.v[(X1.d2+1):(X2.d2+X1.d2)],
                      	     sign(start.v[X2.d2+X1.d2+1])*atanh(i.rho))
                for(i in 1:K) {names(start.v)[i] <- paste("re",i,sep="")
                               names(start.v)[(X1.d2+K+i)] <- paste("re",i,sep="")
                              }	
                nsv <- names(start.v)              
                              
        uid <- unique(id) 
      	N <- length(uid) 
      	uidf <- array(0,N) 
      	for (j in 1:N) uidf[j] <- sum(id==unique(id)[j]) 
      	cuid <- c(0,cumsum(uidf))

      	masses <- rep(1/K,K) 
      	func.opt <- bprobNP

          dat1p <- as.matrix(dat[,3:(X1.d2+2)])
          dat2p <- as.matrix(dat[,(X1.d2+3):(X1.d2+X2.d2+2)])

          dat1 <- function(u) cbind(matrix(replace(rep(0,K),u,1),ncol=K,nrow=n,byrow=TRUE),dat1p)
          dat2 <- function(u) cbind(matrix(replace(rep(0,K),u,1),ncol=K,nrow=n,byrow=TRUE),dat2p)

        if(e.npRE==TRUE){


        T.sv <- func.opt(start.v, dat=dat, dat1=dat1, dat2=dat2, dat1p=dat1p, dat2p=dat2p, 
                         X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, p.weights=p.weights, 
                         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)
                         
          count.npRE <- 0 

      	  while( max(abs(T.sv$gradient)) > pr.tol*10000 && sum(as.numeric(T.sv$gradient=="NaN"))==0){
            
            starm.v <- start.v 
	    start.v <- start.v - ginv(T.sv$hessian)%*%T.sv$gradient 
	    masses <- T.sv$masses
	    T.sv <- func.opt(start.v, dat=dat, dat1=dat1, dat2=dat2, dat1p=dat1p, dat2p=dat2p,  
                             X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, p.weights=p.weights, 
                             K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)
	    count.npRE <- count.npRE + 1
	    
	    if(count.npRE > 5000) break
                                                                           
            }
            if(sum(as.numeric(T.sv$gradient=="NaN"))==0) start.v <- as.vector(start.v) else start.v <- as.vector(starm.v) 
            names(start.v) <- nsv   

                        }                    
                                  
                    } 

                     }else{

  if(npRE==TRUE) stop("Routine for bivariate probit sample selection modelling with random effects not implemented yet")
  inde <- gam1$y > 0
  environment(formula.eq2) <- environment(NULL)
  gam2  <- gam(formula.eq2, binomial(link="probit"), gamma=gamma, weights=p.weights, data=data, subset=inde)
  environment(gam2$formula) <- environment(gam1$formula)
  X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]; X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2)
  y2 <- rep(0,length(inde)); y2[inde] <- gam2$y
  l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth)
  dat <- cbind(gam1$y,y2,X1,X2); n <- length(dat[,1])

  dat1 <- as.matrix(dat[,3:(X1.d2+2)])
  dat2 <- as.matrix(dat[,(X1.d2+3):(X1.d2+X2.d2+2)])

  if(is.null(start.v)) startvSS <- startSS(gam1, gam2, formula.eq2, data, gamma, p.weights, inde, l.sp1, l.sp2, fp); start.v <- startvSS$start.v 
		    	    
		    
  if(H==TRUE && l.sp1==0 && l.sp2==0) func.opt <- bprobSSH else func.opt <- bprobSS		    
  sp <- c(gam1$sp,startvSS$gam2.1$sp)
  

                          }

    if(l.sp1!=0 && l.sp2!=0 && fp==FALSE && npRE==TRUE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2,K=K,npRE)}
    if(l.sp1!=0 && l.sp2!=0 && fp==FALSE && npRE!=TRUE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2,K=K,npRE=FALSE)}
      

    fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, dat=dat, dat1=dat1, dat2=dat2, dat1p=dat1p, dat2p=dat2p, 
                  X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, p.weights=p.weights, 
                  fterm=fterm, mterm=mterm, iterlim = 1e+4, 
                  K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)  

    iter.if <- fit$iterations    
           
    if(fp==TRUE) aut.sp <- FALSE      
           
    if(aut.sp==TRUE){

      if(l.sp1!=0 && l.sp2!=0){
      
       j.it <- stoprule.SP <- 1; conv.sp <- TRUE; count.npREsp <- 0  

	  while( stoprule.SP > pr.tol ){ 

             So <- spS(sp,gam1,gam2); coefo <- fit$argument   
  
		 if(npRE!=TRUE) wor.c <- try(working.comp(fit,X1,X2,X1.d2,X2.d2,n)) else wor.c <- try(working.compNP(fit,X1,X2,X1.d2,X2.d2,n,K,dat1,dat2)) 
                 if(class(wor.c)=="try-error") break
             
                	bs.mgfit <- try(magic(y=wor.c$rW.Z,X=wor.c$rW.X,sp=sp,S=qu.mag$Ss,
                        	           off=qu.mag$off,rank=qu.mag$rank,n.score=3*n,
                                	   gcv=FALSE,gamma=gamma,control=control))
                	if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                	sp <- bs.mgfit$sp
                               
             S <- spS(sp,gam1,gam2)
             o.ests <- c(fit$argument)
            

               if(npRE==TRUE && e.npREsp==TRUE){

                     noe <- names(o.ests)
                     T.sv <- func.opt(o.ests, dat=dat, dat1=dat1, dat2=dat2, dat1p=dat1p, dat2p=dat2p,  
                                      X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, p.weights=p.weights, 
                         	      K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)
                         
      	  		while( max(abs(T.sv$gradient)) > pr.tol*10000 && sum(as.numeric(T.sv$gradient=="NaN"))==0){
      
                        o.estm <- o.ests
	    		o.ests <- o.ests - ginv(T.sv$hessian)%*%T.sv$gradient
	    		masses <- T.sv$masses
	    		T.sv <- func.opt(o.ests, dat=dat, dat1=dat1, dat2=dat2, dat1p=dat1p, dat2p=dat2p,  
            		                 X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, p.weights=p.weights, 
                		         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)
	    		count.npREsp <- count.npREsp + 1
	    		if(count.npREsp > 5000) break
                		                                       }

            		if(sum(as.numeric(T.sv$gradient=="NaN"))==0) o.ests <- as.vector(o.ests) else o.ests <- as.vector(o.estm) 
			names(o.ests) <- noe; fit$argument <- o.ests   
                             }


             fit <- try(trust(func.opt, fit$argument, rinit=rinit, rmax=rmax, dat=dat, dat1=dat1, dat2=dat2, dat1p=dat1p, dat2p=dat2p,  
                              X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, p.weights=p.weights, 
                              iterlim=1e+4, fterm=fterm, mterm=mterm, K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses),silent=TRUE)

              if(class(fit)=="try-error"){ 
               fit  <- trust(func.opt, coefo, rinit=rinit, rmax=rmax, dat=dat, dat1=dat1, dat2=dat2, dat1p=dat1p, dat2p=dat2p, 
                             X1.d2=X1.d2, X2.d2=X2.d2, S=So, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, p.weights=p.weights,
                             iterlim=1e+4, fterm=fterm, mterm=mterm, K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)
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
        T.sv <- func.opt(fit$argument, dat=dat, dat1=dat1, dat2=dat2, dat1p=dat1p, dat2p=dat2p, 
                         X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, p.weights=p.weights,
                         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses)

        npar <- K + X1.d2 + K + X2.d2 + 1 + K - 1
        H.cor <- G.cor <- matrix(0,nrow=npar,ncol=npar)

        for (i in 1:N){
              h <- matrix(0,nrow=npar,ncol=1)
           for (l in 1:K){

              c1 <- T.sv$dl.dbe1[(cuid[i]+1):(cuid[i+1]),l]*dat1(l)[(cuid[i]+1):(cuid[i+1]),]
              if (uidf[i]>1) c1<-apply(c1,2,sum)

              c2 <- T.sv$dl.dbe2[(cuid[i]+1):(cuid[i+1]),l]*dat2(l)[(cuid[i]+1):(cuid[i+1]),]
              if (uidf[i]>1) c2<-apply(c2,2,sum)

              c3 <- sum(T.sv$dl.drho[(cuid[i]+1):(cuid[i+1]),l])

              if (l < K) {c4 <- array(0,K-1); c4[l] <- (1/T.sv$masses[l])}else{ c4 <- array((-1/T.sv$masses[K]),K-1)}

               g <- as.matrix(c(c1,c2,c3,c4))

              H.cor <- H.cor + T.sv$W[i,l]*tcrossprod(g)

              h <- h + T.sv$W[i,l]*g 
           }
           G.cor <- G.cor + tcrossprod(h)
        }

        I.prob1 <- apply(t((1/T.sv$masses^2)*t(T.sv$W)),2,sum)
        if(K!=2) I.prob  <- diag(I.prob1[1:K-1]) + matrix(I.prob1[K],K-1,K-1) else I.prob  <- sum(I.prob1) # I.prob1[1:K-1] + I.prob1[K] 

        He <- adiag(T.sv$hessian,I.prob) - H.cor + G.cor   

        logL <- sum(log(apply(T.sv$Wp3,1,sum))) 

        u1 <- fit$argument[1:K]
        u2 <- fit$argument[(K+X1.d2+1):(K+X1.d2+K)]

        nw <- T.sv$Wp3/matrix(rep(apply(T.sv$Wp3,1,sum),K),ncol=K)

        eb.u1 <- apply(t(u1*t(nw)),1,sum)
        eb.u2 <- apply(t(u2*t(nw)),1,sum)

        Eb.u1 <- rep(eb.u1,uidf)
        Eb.u2 <- rep(eb.u2,uidf)

        fit$eta1 <- Eb.u1 + fit$eta1 
        fit$eta2 <- Eb.u2 + fit$eta2 

    }else{
    
    He <- fit$hessian
    logL <- fit$l}
    
    He.eig <- eigen(He)
    k.e <- sum(as.numeric(He.eig$val<sqrt(.Machine$double.eps)))
    
     if(k.e!=0){
      ind.e <- (length(He.eig$val)-(k.e-1)):length(He.eig$val)
      min.e <- min(He.eig$val[1:(ind.e[1]-1)])
      for(i in 1:k.e) He.eig$val[ind.e[i]] <- min.e/10^i  
      Vb <- He.eig$vec%*%tcrossprod(diag(1/He.eig$val),He.eig$vec)      
     }else{
      Vb <- He.eig$vec%*%tcrossprod(diag(1/He.eig$val),He.eig$vec) 
    }
                                              
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){ if(npRE==TRUE) fit$S.h <- adiag(fit$S.h,matrix(0,K-1,K-1)) 
                                             HeSh <- He - fit$S.h; F <- Vb%*%HeSh
                                       }else{HeSh <- He; F <- Vb%*%HeSh}      
  t.edf <- sum(diag(F))
  rho <- tanh(fit$argument[length(fit$argument)])

  if(selection==TRUE){

  eta1      <- fit$eta1
  non.sel.d <- X1[,names(gam2$coef)]
  param     <- fit$argument[-c(1:length(gam1$coef),length(fit$argument))]
  eta2      <- non.sel.d%*%param

  p00 <- pnorm2(-eta1,-eta2,rho)
  p01 <- pnorm2(-eta1,eta2,-rho)

  fit$p00 <- p00 
  fit$p01 <- p01

  }


L <- list(fit=fit, gam1=gam1, gam2=gam2, gam2.1=startvSS$gam2.1, p.weights=p.weights, sp=sp, iter.sp=j.it, iter.if=iter.if,
          rho=rho, n=n, n.sel=length(gam2$y), X1=X1, X2=X2, X1.d2=X1.d2, X2.d2=X2.d2, 
          l.sp1=l.sp1, l.sp2=l.sp2, He=He, HeSh=HeSh, Vb=Vb, F=F, 
          t.edf=t.edf, bs.mgfit=bs.mgfit, conv.sp=conv.sp, wor.c=wor.c,
          p11=fit$p11, p10=fit$p10, p01=fit$p01, p00=fit$p00, eta1=fit$eta1, eta2=fit$eta2, 
          dat=dat,sel=selection,masses=fit$masses,K=K,iter.npRE=count.npRE, iter.npREsp=count.npREsp, 
          npRE=npRE, logL=logL, eb.u1=eb.u1, eb.u2=eb.u2, Eb.u1=Eb.u1, Eb.u2=Eb.u2, id=id, uidf=uidf, T.sv=T.sv)

class(L) <- "SemiParBIVProbit"

L

}

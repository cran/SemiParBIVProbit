startSS <- function(gam1, gam2, formula.eq2, data, gamma, inde, l.sp1, l.sp2, fp){

                               p.g1 <- predict(gam1)
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

list(start.v=start.v,gam2.1=gam2.1)

}


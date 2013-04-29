startSS <- function(gam1, gam2, formula.eq2, data, gamma, weights, inde, l.sp1, l.sp2){

    p.g1 <- predict(gam1)
    imr <- data$imr <- dnorm(p.g1)/pnorm(p.g1)
    formula.eq2.1 <- update.formula(formula.eq2, ~. + imr)
    ols <- eval(substitute(gam(formula.eq2.1, data=data, gamma=gamma, weights=weights, subset=inde),list(weights=weights,inde=inde)))
    rhi <- coef(ols)["imr"]/sqrt(ols$sig2)
    ta2 <- (1 + rhi^2*imr*(-p.g1-imr))[inde]

    M <- ols$model; vSn <- names(M) 
    fw <-  paste(vSn[1],"~",vSn[2],sep="") 
    if(dim(M)[2]>2) for (i in 3:(length(vSn)-1)) fw <- paste(fw,"+",vSn[i],sep="")
    fw <- as.formula(fw)  

    MM <- as.data.frame(model.matrix(fw,M)[,-1])
                               
    if(dim(M)[2]==2) names(MM) <- vSn[2]
    data.r <- MM/sqrt( pmax(10000*.Machine$double.eps, ta2) ) 
    data.r <- as.data.frame(cbind(ols$y,data.r)); names(data.r)[1] <- names(M)[1]
                        
    if(l.sp2!=0) l <- as.matrix(data.r[,2:(ols$smooth[[1]]$first.para-1)]) else l <- as.matrix(data.r[,-1])

    fw <-  paste(names(data.r)[1],"~ l",sep="")

    if(l.sp2!=0) for(i in 1:l.sp2) fw <- paste(fw,"+",ols$smooth[[i]]$label,sep="")

    fw <- as.formula(fw)  

gam2.1 <- suppressWarnings(try(eval(substitute(gam(fw, binomial(link="probit"), gamma=gamma, weights=weights[inde], data=data.r),list(weights=weights,inde=inde))),silent=TRUE))

                               environment(gam2.1$formula) <- environment(gam2$formula) 

                               if(class(gam2.1)[1]=="try-error"){
                                                                 names(rhi) <- "athrho" 
                                                                 start.v <- c(coef(gam1),coef(gam2),atanh(rhi))
                                                                 gam2.1 <- gam2 
                                                                } else {

                               if(l.sp2!=0) names(gam2.1$coefficients)[1:(ols$smooth[[1]]$first.para-1)] <- names(ols$coefficients)[1:(ols$smooth[[1]]$first.para-1)] 
                               else names(gam2.1$coefficients) <- names(ols$coefficients) 

                               rho.c <- gam2.1$coefficients["imr"]  
                               rho <- ifelse( abs(rho.c) > 0.90, sign(rho.c)*0.90, rho.c ); names(rho) <- "athrho"
  			       start.v <- c(coef(gam1),gam2.1$coef[names(gam2.1$coefficients)!="imr"], atanh(rho)  )
                                                                       }
list(start.v=start.v,gam2.1=gam2.1)

}



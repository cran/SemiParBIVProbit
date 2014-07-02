summary.SemiParBIVProbit <- function(object,n.sim=1000,s.meth="svd",prob.lev=0.05,thrs1=0.5,thrs2=0.5,...){

  testStat <- function (p, X, V, rank = NULL) {
      qrx <- qr(X)
      R <- qr.R(qrx)
      V <- R %*% tcrossprod(V[qrx$pivot, qrx$pivot],R)
      V <- (V + t(V))/2
      ed <- eigen(V, symmetric = TRUE)
      k <- max(0, floor(rank))
      nu <- abs(rank - k)

          if (rank > k + 0.05 || k == 0) 
              k <- k + 1
          nu <- 0
          rank <- k
      
      if (nu > 0) 
          k1 <- k + 1
      else k1 <- k
      r.est <- sum(ed$values > max(ed$values) * .Machine$double.eps^0.9)
      if (r.est < k1) {
          k1 <- k <- r.est
          nu <- 0
          rank <- r.est
      }
      vec <- ed$vectors
      if (k1 < ncol(vec)) 
          vec <- vec[, 1:k1, drop = FALSE]
      if (k == 0) {
          vec <- t(t(vec) * sqrt(nu/ed$val[1]))
      }
      if (nu > 0 && k > 0) {
          if (k > 1) 
              vec[, 1:(k - 1)] <- t(t(vec[, 1:(k - 1)])/sqrt(ed$val[1:(k - 
                  1)]))
          b12 <- 0.5 * nu * (1 - nu)
          if (b12 < 0) 
              b12 <- 0
          b12 <- sqrt(b12)
          B <- matrix(c(1, b12, b12, nu), 2, 2)
          ev <- diag(ed$values[k:k1]^-0.5)
          B <- ev %*% B %*% ev
          eb <- eigen(B, symmetric = TRUE)
          rB <- eb$vectors %*% tcrossprod(diag(sqrt(eb$values)),eb$vectors)
          vec[, k:k1] <- t(tcrossprod(rB,vec[, k:k1]))
      }
      else {
          vec <- t(t(vec)/sqrt(ed$val[1:k]))
      }
      d <- crossprod(vec,R%*%p)
      d <- sum(d^2)
      attr(d, "rank") <- rank
      d
}

  good <- object$fit$good
  n.good <- sum(as.numeric(good))
  tableN <- list(NULL,NULL)
  table <- list()
  n.sel <- object$n.sel; masses <- object$masses
  CIl1 <- CIl2 <- table.R <- table.P <- table.F <- P1 <- P2 <- QPS1 <- QPS2 <- CR1 <- CR1 <- CR2 <- MR <- CIkt <- CId <- NULL  
  epsilon <- .Machine$double.eps*10^6
  est.RHOb <- est.KeTb <- est.l1 <- est.l2 <- ass.ps2 <- rep(NA,n.sim) 
  bb <- c("BB1.0","BB1.180","BB1.90","BB1.270",
          "BB6.0","BB6.180","BB6.90","BB6.270",
          "BB7.0","BB7.180","BB7.90","BB7.270",
          "BB8.0","BB8.180","BB8.90","BB8.270")
 
  lf.n <- lf <- length(coef(object))
  F  <- object$F[1:lf,1:lf]
  Vr <- object$Vb[1:lf,1:lf] 

  #object$fit$hessian

  #Hll=HH[(d-1):d,(d-1):d]
  #Hbl=HH[(d-1):d,1:(d-2)]
  #Hlb=HH[1:(d-2),(d-1):d]
  #Hbb=HH[1:(d-2),1:(d-2)]
  #Hbb-Hlb%*%solve(Hll)%*%Hbl


          
  SE <- sqrt(diag(object$Vb[1:lf,1:lf]))
  n  <- object$n 


  bs <- rmvnorm(n.sim, mean = coef(object), sigma=object$Vb, method=s.meth)

  if(object$BivD %in% bb ) lf.n <- lf-1

  if(object$PL != "P" && object$fitPL!="fixed") { if(object$eqPL=="both") lf.n <- lf-2 else lf.n <- lf-1 }

   if(object$BivD %in% c("N","T"))      {est.RHOb <- tanh(bs[,lf.n ]); est.RHOb <- ifelse(est.RHOb %in% c(-1,1), sign(est.RHOb)*0.9999999, est.RHOb)}
   if(object$BivD=="F")                  est.RHOb <- bs[,lf.n ] + epsilon

   if(object$BivD %in% c("C0", "C180") ) est.RHOb <-   exp(bs[,lf.n ]) + epsilon  
   if(object$BivD %in% c("C90","C270") ) est.RHOb <- -(exp(bs[,lf.n ]) + epsilon) 

   if(object$BivD %in% c("J0", "J180") ) est.RHOb <-   1+exp(bs[,lf.n]) + epsilon  
   if(object$BivD %in% c("J90","J270") ) est.RHOb <- -(1+exp(bs[,lf.n]) + epsilon) 

   if(object$BivD %in% c("G0", "G180") ) est.RHOb <-   1+exp(bs[,lf.n])  
   if(object$BivD %in% c("G90","G270") ) est.RHOb <- -(1+exp(bs[,lf.n])) 
   
   if(object$BivD %in% c("BB1.0", "BB1.180")){est.RHOb <-   exp(bs[,lf.n]) + epsilon;  ass.ps2 <-   exp(bs[,lf]) + 1}
   if(object$BivD %in% c("BB1.90","BB1.270")){est.RHOb <- -(exp(bs[,lf.n]) + epsilon); ass.ps2 <- -(exp(bs[,lf]) + 1)}
   
   if(object$BivD %in% c("BB6.0", "BB6.180")){est.RHOb <-   exp(bs[,lf.n]) + 1; ass.ps2 <-   exp(bs[,lf]) + 1}
   if(object$BivD %in% c("BB6.90","BB6.270")){est.RHOb <- -(exp(bs[,lf.n]) + 1);ass.ps2 <- -(exp(bs[,lf]) + 1)}
   
   if(object$BivD %in% c("BB7.0", "BB7.180")){est.RHOb <-   exp(bs[,lf.n]) + 1; ass.ps2 <-   exp(bs[,lf]) + epsilon}
   if(object$BivD %in% c("BB7.90","BB7.270")){est.RHOb <- -(exp(bs[,lf.n]) + 1);ass.ps2 <- -(exp(bs[,lf]) + epsilon)}
   
   if(object$BivD %in% c("BB8.0", "BB8.180")){est.RHOb <-   exp(bs[,lf.n]) + 1; ass.ps2 <-  pnorm(bs[,lf]); ass.ps2 <- ifelse(ass.ps2==0, epsilon, ass.ps2)}
   if(object$BivD %in% c("BB8.90","BB8.270")){est.RHOb <- -(exp(bs[,lf.n]) + 1);ass.ps2 <- -pnorm(bs[,lf]); ass.ps2 <- ifelse(ass.ps2==0, -epsilon, ass.ps2)}  
   
   
   if(object$BivD=="T") asp2 <- rep(object$nu,n.sim) else asp2 <- ass.ps2
   
   for(i in 1:n.sim) est.KeTb[i] <- BiCopPar2Tau(object$nC,est.RHOb[i],par2=asp2[i]) # this is not efficient but the function gives problems with vectors...

   if((object$PL=="PP" || object$PL=="RPP") && object$fitPL!="fixed"){ if(object$eqPL=="both"){   est.l1 <- exp(bs[, lf.n+1 ]) + epsilon; est.l2 <- exp(bs[,lf ]) + epsilon}
                       if(object$eqPL=="first"){  est.l1 <- exp(bs[, lf ]) + epsilon;     est.l2 <- 1}
                       if(object$eqPL=="second"){ est.l2 <- exp(bs[, lf ]) + epsilon;     est.l1 <- 1}
                     }  
                     
   if(object$PL=="SN" && object$fitPL!="fixed"){ if(object$eqPL=="both"){   est.l1 <- bs[, lf.n+1 ]; est.l2 <- bs[,lf ]}
                       if(object$eqPL=="first"){  est.l1 <- bs[, lf ];     est.l2 <- 0}
                       if(object$eqPL=="second"){ est.l2 <- bs[, lf ];     est.l1 <- 0}
                     }                       

             
  CIrs <- as.numeric(quantile(est.RHOb,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))
  CIkt <- as.numeric(quantile(est.KeTb,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))
  
  if(object$PL != "P" && object$fitPL!="fixed"){
                       CIl1 <- as.numeric(quantile(est.l1,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))
                       CIl2 <- as.numeric(quantile(est.l2,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))
                       }
  if(object$BivD %in% bb ) CId <- as.numeric(quantile(ass.ps2,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))

                     
  ind <- list(ind1=1:(object$gam1$nsdf),ind2=object$X1.d2+(1:(object$gam2$nsdf)))

  for(i in 1:2){
  estimate <- coef(object)[ind[[i]]]
  se       <- SE[ind[[i]]]
  ratio    <- estimate/se
  pv       <- 2*pnorm(abs(ratio), lower.tail = FALSE)
  table[[i]] <- cbind(estimate,se,ratio,pv)
  dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }


  l.sp11 <- length(object$gam1$smooth)
  l.sp22 <- length(object$gam2$smooth) 
  if( (l.sp11!=0 || l.sp22!=0) ){

  	pTerms.df <- pTerms.chi.sq <- pTerms.pv <- edf <- tableN <- list(0,0)
        XX <- cbind(object$X1,object$X2)
        
           for(i in 1:2){

             if(i==1) {mm <- l.sp11; if(mm==0) next}
             if(i==2) {mm <- l.sp22; if(mm==0) break} 
  
		for(k in 1:mm){

                        if(i==1){gam <- object$gam1; ind <- (gam$smooth[[k]]$first.para):(gam$smooth[[k]]$last.para)} 
                        else{gam <- object$gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para)+object$X1.d2} 
			edf[[i]][k] <- sum(diag(F)[ind])
			names(edf[[i]])[k] <- gam$smooth[[k]]$label 
			b  <- coef(object)[ind]
			V  <- Vr[ind,ind]
			Xt <- XX[, ind] 
			pTerms.df[[i]][k] <- min(ncol(Xt), edf[[i]][k])
			pTerms.chi.sq[[i]][k] <- Tp <- testStat(b, Xt, V, pTerms.df[[i]][k])
			pTerms.df[[i]][k] <- attr(Tp, "rank")
                        pTerms.pv[[i]][k] <- pchisq(pTerms.chi.sq[[i]][k], df = pTerms.df[[i]][k], lower.tail = FALSE)
			                 
                }
              tableN[[i]] <- cbind(edf[[i]], pTerms.df[[i]], pTerms.chi.sq[[i]], pTerms.pv[[i]])
              dimnames(tableN[[i]])[[2]] <- c("edf", "Est.rank", "Chi.sq", "p-value")
            }

  }




 if(object$sel==FALSE){
 
 Pre.p <- matrix(NA,n.good,8)
 Pre.c <- matrix(NA,n.good,2)


 Pre.p[,1:6] <- cbind(object$y1[good],object$y2[good],object$p11,object$p10,object$p01,object$p00)

 for(i in 1:n.good) {
   ind <- sort(Pre.p[i,3:6],index.return=TRUE)$ix[4]
   if(ind==1) Pre.p[i,7:8] <- c(1,1) 
   if(ind==2) Pre.p[i,7:8] <- c(1,0) 
   if(ind==3) Pre.p[i,7:8] <- c(0,1) 
   if(ind==4) Pre.p[i,7:8] <- c(0,0) 
 }

 Pre.p <- Pre.p[,-c(3:6)]

 for(i in 1:n.good){
   Pre.c[i,1] <- paste(as.character(Pre.p[i,1:2]), collapse="")
   Pre.c[i,2] <- paste(as.character(Pre.p[i,3:4]), collapse="")
 }

 matches <- as.numeric(Pre.c[,1]==Pre.c[,2])
 MR <- mean(matches)*100

 c00.p <- c10.p <- c01.p <- c11.p <- 0

 for(i in 1:n.good){
  if(Pre.c[i,1]=="00" & Pre.c[i,1]==Pre.c[i,2]) c00.p <- c00.p + 1
  if(Pre.c[i,1]=="10" & Pre.c[i,1]==Pre.c[i,2]) c10.p <- c10.p + 1
  if(Pre.c[i,1]=="01" & Pre.c[i,1]==Pre.c[i,2]) c01.p <- c01.p + 1
  if(Pre.c[i,1]=="11" & Pre.c[i,1]==Pre.c[i,2]) c11.p <- c11.p + 1
 }

 table.R <- as.data.frame(matrix(table(Pre.c[,1]),2,2,byrow=TRUE))
 table.P <- as.data.frame(matrix(c(c00.p,c01.p,c10.p,c11.p),2,2,byrow=TRUE))
 table.F <- table.P/table.R
 dimnames(table.R)[[1]] <- dimnames(table.R)[[2]] <- dimnames(table.P)[[1]] <- dimnames(table.P)[[2]] <- dimnames(table.F)[[1]] <- dimnames(table.F)[[2]] <- c("0", "1")
 
 P1 <- object$p11 + object$p10
 P2 <- object$p11 + object$p01
 QPS1 <- 1/n*( sum( 2*( object$y1[good] - P1)^2 ) )
 QPS2 <- 1/n*( sum( 2*( object$y2[good] - P2)^2 ) )

 P1.b <- ifelse(P1 > thrs1, 1, 0)
 P2.b <- ifelse(P2 > thrs2, 1, 0)

 CR1 <- mean(as.numeric(object$y1[good]==P1.b))*100
 CR2 <- mean(as.numeric(object$y2[good]==P2.b))*100

 }

  
  res <- list(tableP1=table[[1]], tableP2=table[[2]], 
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], 
              n=n, rho=object$rho, theta=object$theta, delta=object$delta, KeT=object$KeT,   
              formula1=object$gam1$formula, formula2=object$gam2$formula, 
              l.sc1=l.sp11, l.sc2=l.sp22, pPen1=object$pPen1, pPen2=object$pPen2,  
              t.edf=object$t.edf, CIrs=CIrs, CIkt=CIkt, CIl1=CIl1, CIl2=CIl2, 
              CId=CId, 
              sel=object$sel,n.sel=n.sel, 
              BivD=object$BivD,nu=object$nu, 
              PL=object$PL, xi1=object$xi1, xi2=object$xi2,
              table.R=table.R, table.P=table.P, table.F=table.F, MR=MR,
              P1=P1, P2=P2, QPS1=QPS1, QPS2=QPS2, CR1=CR1, CR2=CR2,
              good=good
              )
  class(res) <- "summary.SemiParBIVProbit"
      
                                        

res

}




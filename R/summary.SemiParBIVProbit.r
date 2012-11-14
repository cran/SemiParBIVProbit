summary.SemiParBIVProbit <- function(object,n.sim=1000,s.meth="svd",sig.lev=0.05,thrs1=0.5,thrs2=0.5,...){

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
  
  tableN <- list(NULL,NULL)
  n.sel <- object$n.sel; masses <- object$masses
  table.R <- table.P <- table.F <- P1 <- P2 <- QPS1 <- QPS2 <- CR1 <- CR1 <- CR2 <- MR <- table.npRE <- NULL  

  lf <- length(object$fit$argument)
  F  <- object$F[1:lf,1:lf]
  Vr <- object$Vb[1:lf,1:lf] 
          
  SE <- sqrt(diag(object$Vb[1:lf,1:lf]))
  n  <- object$n 

  if(object$npRE==FALSE) bs <- rmvnorm(n.sim, mean = object$fit$argument, sigma=object$Vb, method=s.meth)
  else bs <- rmvnorm(n.sim, mean = c(object$fit$argument,object$fit$masses[1:(object$K-1)]), sigma=object$Vb, method=s.meth)

  est.RHOb <- rep(NA,n.sim)
  for(i in 1:n.sim) est.RHOb[i] <- tanh(bs[i,lf])
  CIrs <- as.numeric(quantile(est.RHOb,c(sig.lev/2,1-sig.lev/2)))
  
  Kk1 <- Kk2 <- 0; if(object$npRE==TRUE){ Kk1 <- object$K - 1; Kk2 <- Kk1 + 1}  

  table <- list()
  ind <- list(ind1=1:(object$gam1$nsdf+Kk1),ind2=object$X1.d2+Kk2+(1:(object$gam2$nsdf+Kk1)))

  for(i in 1:2){
  estimate <- object$fit$argument[ind[[1]]]
  se       <- SE[ind[[1]]]
  ratio    <- estimate/se
  pv       <- 2*pnorm(abs(ratio), lower.tail = FALSE)
  table[[i]] <- cbind(estimate,se,ratio,pv)
  dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }

  if(object$npRE==TRUE){
  	table.npRE <- cbind(table[[1]][c(1:Kk2),1],table[[2]][c(1:Kk2),1],object$masses)
  	dimnames(table.npRE)[[2]] <- c("Eq. 1", "Eq. 2", "masses")
  	table[[1]] <- table[[1]][-c(1:Kk2),]; table[[2]] <- table[[2]][-c(1:Kk2),] 
                       }

  if(object$l.sp1!=0 && object$l.sp2!=0){

  	pTerms.df <- pTerms.chi.sq <- pTerms.pv <- edf <- tableN <- list(0,0)
        
           for(i in 1:2){
                if(i==1) mm <- object$l.sp1 else mm <- object$l.sp2
  
		for(k in 1:mm){

                        if (i==1) {gam <- object$gam1; ind <- (object$gam1$smooth[[k]]$first.para+Kk1):(object$gam1$smooth[[k]]$last.para+Kk1)} else{gam <- object$gam2; ind <- (object$gam2$smooth[[k]]$first.para:object$gam2$smooth[[k]]$last.para+Kk1)+object$X1.d2+Kk2}
			edf[[i]][k] <- sum(diag(F)[ind])
			names(edf[[i]])[k] <- gam$smooth[[k]]$label 
			b  <- object$fit$argument[ind]
			V  <- Vr[ind,ind]
			if(i==1) Xt <- object$X1[, 1:length(ind)+gam$nsdf] else Xt <- object$X2[, 1:length(ind)+gam$nsdf]
			pTerms.df[[i]][k] <- min(ncol(Xt), edf[[i]][k])
			pTerms.chi.sq[[i]][k] <- Tp <- testStat(b, Xt, V, pTerms.df[[i]][k])
			pTerms.df[[i]][k] <- attr(Tp, "rank")
                        pTerms.pv[[i]][k] <- pchisq(pTerms.chi.sq[[i]][k], df = pTerms.df[[i]][k], lower.tail = FALSE)
			                 
                }
              tableN[[i]] <- cbind(edf[[i]], pTerms.df[[i]], pTerms.chi.sq[[i]], pTerms.pv[[i]])
              dimnames(tableN[[i]])[[2]] <- c("edf", "Est.rank", "Chi.sq", "p-value")
            }

  }


 if(object$npRE==FALSE){ 

 if(object$sel==FALSE){
 
 Pre.p <- matrix(NA,n,8)
 Pre.c <- matrix(NA,n,2)

 Pre.p[,1:6] <- cbind(object$dat[,1:2],object$p11,object$p10,object$p01,object$p00)

 for(i in 1:n) {
   ind <- sort(Pre.p[i,3:6],index.return=TRUE)$ix[4]
   if(ind==1) Pre.p[i,7:8] <- c(1,1) 
   if(ind==2) Pre.p[i,7:8] <- c(1,0) 
   if(ind==3) Pre.p[i,7:8] <- c(0,1) 
   if(ind==4) Pre.p[i,7:8] <- c(0,0) 
 }

 Pre.p <- Pre.p[,-c(3:6)]

 for(i in 1:n){
   Pre.c[i,1] <- paste(as.character(Pre.p[i,1:2]), collapse="")
   Pre.c[i,2] <- paste(as.character(Pre.p[i,3:4]), collapse="")
 }

 matches <- as.numeric(Pre.c[,1]==Pre.c[,2])
 MR <- mean(matches)*100

 c00.p <- c10.p <- c01.p <- c11.p <- 0

 for(i in 1:n){
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
 QPS1 <- 1/n*( sum( 2*( object$gam1$y - P1)^2 ) )
 QPS2 <- 1/n*( sum( 2*( object$gam2$y - P2)^2 ) )

 P1.b <- ifelse(P1 > thrs1, 1, 0)
 P2.b <- ifelse(P2 > thrs2, 1, 0)

 CR1 <- mean(as.numeric(object$gam1$y==P1.b))*100
 CR2 <- mean(as.numeric(object$gam2$y==P2.b))*100

 }
 }

  
  res <- list(tableP1=table[[1]], tableP2=table[[2]], 
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], 
              n=n, rho=object$rho, 
              formula1=object$gam1$formula, formula2=object$gam2$formula, 
              l.sp1=object$l.sp1, l.sp2=object$l.sp2, 
              t.edf=object$t.edf, CIrs=CIrs, sel=object$sel,n.sel=n.sel, npRE=object$npRE,
              table.R=table.R, table.P=table.P, table.F=table.F, MR=MR,
              P1=P1, P2=P2, QPS1=QPS1, QPS2=QPS2, CR1=CR1, CR2=CR2,
              masses=object$masses, table.npRE=table.npRE
              )
  class(res) <- "summary.SemiParBIVProbit"
      
                                        

res

}




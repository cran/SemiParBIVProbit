summary.SemiParBIVProbit <- function(object,n.sim=1000,s.meth="svd",sig.lev=0.05,thrs1=0.5,thrs2=0.5,...){

  testStat <- function (p, X, V, rank = NULL) {
      qrx <- qr(X)
      R <- qr.R(qrx)
      V <- R %*% V[qrx$pivot, qrx$pivot] %*% t(R)
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
          rB <- eb$vectors %*% diag(sqrt(eb$values)) %*% t(eb$vectors)
          vec[, k:k1] <- t(rB %*% t(vec[, k:k1]))
      }
      else {
          vec <- t(t(vec)/sqrt(ed$val[1:k]))
      }
      d <- t(vec) %*% (R %*% p)
      d <- sum(d^2)
      attr(d, "rank") <- rank
      d
}
  
  lf <- length(object$fit$argument)
  F  <- object$F[1:lf,1:lf]
  Vr <- object$Vb[1:lf,1:lf] 
          
  SE <- sqrt(diag(object$Vb[1:lf,1:lf]))
  n  <- object$n 

  #if(object$npRE==TRUE){ mm <- c(object$fit$argument,object$fit$masses[1:2])
  #    bs <- rmvnorm(n.sim, mean = mm, sigma=object$Vb, method=s.meth)}else{
  bs <- rmvnorm(n.sim, mean = object$fit$argument, sigma=object$Vb[1:lf,1:lf], method=s.meth)#}
  d.rho <- dim(object$Vb[1:lf,1:lf])[1]
  est.RHOb <- rep(NA,n.sim)
  for(i in 1:n.sim) est.RHOb[i] <- tanh(bs[i,d.rho])
  CIrs <- as.numeric(quantile(est.RHOb,c(sig.lev/2,1-sig.lev/2)))
  
  Kk1 <- Kk2 <- 0; if(object$npRE==TRUE){ Kk1 <- object$K - 1; Kk2 <- Kk1 + 1}  

  estimate1 <- object$fit$argument[1:(object$gam1$nsdf+Kk1)]
  se1       <- SE[1:(object$gam1$nsdf+Kk1)]
  ratio1    <- estimate1/se1
  pv1       <- 2*pnorm(abs(ratio1), lower.tail = FALSE)
  table1    <- cbind(estimate1,se1,ratio1,pv1)

  estimate2 <- object$fit$argument[object$X1.d2+Kk2+(1:(object$gam2$nsdf+Kk1))]
  se2       <- SE[object$X1.d2+Kk2+(1:(object$gam2$nsdf+Kk1))]
  ratio2    <- estimate2/se2
  pv2       <- 2*pnorm(abs(ratio2), lower.tail = FALSE)
  table2    <- cbind(estimate2,se2,ratio2,pv2)

  dimnames(table1)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  dimnames(table2)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  if(object$npRE==TRUE){
  	table.npRE <- cbind(table1[c(1:Kk2),1],table2[c(1:Kk2),1],object$masses)
  	dimnames(table.npRE)[[2]] <- c("Eq. 1", "Eq. 2", "masses")
  	table1 <- table1[-c(1:Kk2),]; table2 <- table2[-c(1:Kk2),] 
                       }

  if(object$l.sp1!=0 && object$l.sp2!=0){
  	pTerms.df1 <- pTerms.chi.sq1 <- pTerms.pv1 <- edf1 <- NA
  	pTerms.df2 <- pTerms.chi.sq2 <- pTerms.pv2 <- edf2 <- NA
		for(k in 1:object$l.sp1){
			ind <- (object$gam1$smooth[[k]]$first.para+Kk1):(object$gam1$smooth[[k]]$last.para+Kk1)
			edf1[k] <- sum(diag(F)[ind])
			names(edf1)[k] <- object$gam1$smooth[[k]]$label 
			b  <- object$fit$argument[ind]
			V  <- Vr[ind,ind]
			Xt <- object$X1[, 1:length(ind)+object$gam1$nsdf]
			pTerms.df1[k] <- min(ncol(Xt), edf1[k])
			pTerms.chi.sq1[k] <- Tp <- testStat(b, Xt, V, pTerms.df1[k])
			pTerms.df1[k] <- attr(Tp, "rank")
                        pTerms.pv1[k] <- pchisq(pTerms.chi.sq1[k], df = pTerms.df1[k], lower.tail = FALSE)
			                 
            }
  	table1.1 <- cbind(edf1, pTerms.df1, pTerms.chi.sq1, pTerms.pv1)

		for(k in 1:object$l.sp2){
			ind <- (object$gam2$smooth[[k]]$first.para:object$gam2$smooth[[k]]$last.para+Kk1)+object$X1.d2+Kk2
			edf2[k] <- sum(diag(F)[ind])
			names(edf2)[k] <- object$gam2$smooth[[k]]$label 
			b  <- object$fit$argument[ind]
			V  <- Vr[ind,ind]
			Xt <- object$X2[, 1:length(ind)+object$gam2$nsdf]
			pTerms.df2[k] <- min(ncol(Xt), edf2[k])
			pTerms.chi.sq2[k] <- Tp <- testStat(b, Xt, V, pTerms.df2[k])   
			pTerms.df2[k] <- attr(Tp, "rank")
			pTerms.pv2[k] <- pchisq(pTerms.chi.sq2[k], df = pTerms.df2[k], lower.tail = FALSE)
            }
	table2.2 <- cbind(edf2, pTerms.df2, pTerms.chi.sq2, pTerms.pv2)
  dimnames(table1.1)[[2]] <- c("edf", "Est.rank", "Chi.sq", "p-value")
  dimnames(table2.2)[[2]] <- c("edf", "Est.rank", "Chi.sq", "p-value")
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

  

 if(object$sel==FALSE && object$npRE==FALSE){


  if(object$l.sp1!=0 && object$l.sp2!=0){res <- list(tableP1=table1, tableP2=table2, 
                                           tableNP1=table1.1, tableNP2=table2.2, 
                                           n=n, rho=object$rho, 
                                           formula1=object$gam1$formula, formula2=object$gam2$formula, 
                                           l.sp1=object$l.sp1, l.sp2=object$l.sp2, 
                                           t.edf=object$t.edf, CIrs=CIrs, sel=object$sel,npRE=object$npRE,
                   			   table.R=table.R, table.P=table.P, table.F=table.F, MR=mean(matches)*100,
                   			   P1=P1, P2=P2, QPS1=QPS1, QPS2=QPS2, CR1=CR1, CR2=CR2
                                           )
                               class(res) <- "summary.SemiParBIVProbit"
  }else{res <- list(tableP1=table1,tableP2=table2, 
                   n=n, rho=object$rho, 
                   formula1=object$gam1$formula, formula2=object$gam2$formula, 
                   l.sp1=0, l.sp2=0, npRE=object$npRE,
                   t.edf=object$t.edf, CIrs=CIrs, sel=object$sel,
                   table.R=table.R, table.P=table.P, table.F=table.F, MR=mean(matches)*100,
                   P1=P1, P2=P2, QPS1=QPS1, QPS2=QPS2, CR1=CR1, CR2=CR2
                   )
       class(res) <- "summary.SemiParBIVProbit"
         }
         
                                        }

 if(object$sel==TRUE && object$npRE==FALSE){


 if(object$l.sp1!=0 && object$l.sp2!=0){res <- list(tableP1=table1, tableP2=table2, 
                                           tableNP1=table1.1, tableNP2=table2.2, 
                                           n=n, rho=object$rho, 
                                           formula1=object$gam1$formula, formula2=object$gam2$formula, 
                                           l.sp1=object$l.sp1, l.sp2=object$l.sp2, npRE=object$npRE,
                                           t.edf=object$t.edf, CIrs=CIrs, sel=object$sel, n.sel=object$n.sel, npRE=object$npRE
                                           )
                               class(res) <- "summary.SemiParBIVProbit"
  }else{res <- list(tableP1=table1,tableP2=table2, 
                   n=n, rho=object$rho, 
                   formula1=object$gam1$formula, formula2=object$gam2$formula, 
                   l.sp1=0, l.sp2=0,npRE=object$npRE, 
                   t.edf=object$t.edf, CIrs=CIrs, sel=object$sel, npRE=object$npRE, n.sel=object$n.sel
                   )
       class(res) <- "summary.SemiParBIVProbit"
       }
}



 if(object$npRE==TRUE){


 if(object$l.sp1!=0 && object$l.sp2!=0){res <- list(tableP1=table1, tableP2=table2, 
                                           tableNP1=table1.1, tableNP2=table2.2, 
                                           n=n, rho=object$rho, 
                                           formula1=object$gam1$formula, formula2=object$gam2$formula, 
                                           l.sp1=object$l.sp1, l.sp2=object$l.sp2, 
                                           t.edf=object$t.edf, CIrs=CIrs, sel=object$sel, 
                                           npRE=object$npRE, masses=object$masses, table.npRE=table.npRE
                                           )
                               class(res) <- "summary.SemiParBIVProbit"
  }else{res <- list(tableP1=table1,tableP2=table2, 
                   n=n, rho=object$rho, 
                   formula1=object$gam1$formula, formula2=object$gam2$formula, 
                   l.sp1=0, l.sp2=0, 
                   t.edf=object$t.edf, CIrs=CIrs, sel=object$sel, 
                   npRE=object$npRE, masses=object$masses, table.npRE=table.npRE
                   )
       class(res) <- "summary.SemiParBIVProbit"
       }

                    }

res

}




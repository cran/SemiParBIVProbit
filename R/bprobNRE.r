bprobNRE <- function(params, y1, y2, q1, q2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, weights=weights, X1.d2, X2.d2, sp=NULL, qu.mag=NULL, gp1, gp2, fp, 
                      l.sp1, l.sp2, K=NULL, n, N, cuid, uidf, masses=NULL, NGQ, dat1all, dat2all, W){

  sigma1.st <- params[X1.d2+1]
  corr.st.u <- params[X1.d2+X2.d2+2]
  sigma2.st <- params[X1.d2+X2.d2+3]
  corr.st.e <- params[X1.d2+X2.d2+4]

  sigma1 <- exp(sigma1.st)
  corr.u <- tanh(corr.st.u)
  sigma2 <- exp(sigma2.st)
  corr.e <- tanh(corr.st.e)

  eta1 <- dat1all%*%c(params[1:X1.d2],sigma1)
  eta2 <- dat2all%*%c(params[(X1.d2+2):(X1.d2+X2.d2+1)],corr.u*sigma2,sqrt(pmax(10000*.Machine$double.eps, 1-corr.u^2))*sigma2)
  
## Other definitions (Greene p. 738-741)

  corr.sq <- q1*q2*corr.e
  w1 <- q1*eta1
  w2 <- q2*eta2
  delta <- 1/sqrt(pmax(10000*.Machine$double.eps,1-corr.sq^2))
  v1 <- delta*(w2-corr.sq*w1)
  v2 <- delta*(w1-corr.sq*w2)
  g1 <- dnorm(w1)*pnorm(v1)
  g2 <- dnorm(w2)*pnorm(v2)
  phi2 <- dnorm2(w1,w2,rho=corr.sq)
  
## Leading to log-likelihood (l.par) and normalized weights (WN2)

  PHI2 <- pmax(abs(pnorm2( w1,w2,cov12=corr.sq)),1000*.Machine$double.eps)
  l.par1 <- weights*log(PHI2)

  WN <- matrix(0,nrow=N,ncol=NGQ^2)
  for (r in 1:N)
      for (c in 1:NGQ^2)
          WN[r,c] <- exp(sum(l.par1[seq(c+cuid[r]*NGQ^2,cuid[r+1]*NGQ^2,NGQ^2)]))*W[c]

  l.par <- sum(log(apply(WN,1,sum)))

  WN2 <- (1/apply(WN,1,sum))*WN

## Derivatives for chain rule

  drh.drh.st.e <- 4*exp(2*corr.st.e)/(exp(2*corr.st.e)+1)^2
  drh.drh.st.u <- 4*exp(2*corr.st.u)/(exp(2*corr.st.u)+1)^2
  drh.drh.st2.e <- 8*exp(2*corr.st.e)*(1-exp(2*corr.st.e))/(exp(2*corr.st.e)+1)^3
  drh.drh.st2.u <- 8*exp(2*corr.st.u)*(1-exp(2*corr.st.u))/(exp(2*corr.st.u)+1)^3
  dsig1.dsig1.st <- exp(sigma1.st)
  dsig2.dsig2.st <- exp(sigma2.st)
  
  

  TbT1 <- sigma2*drh.drh.st.u
  TbT2 <- corr.u*dsig2.dsig2.st
  TbT3 <- -sigma2*corr.u*drh.drh.st.u/sqrt(pmax(10000*.Machine$double.eps, 1-corr.u^2))
  TbT4 <- sqrt(pmax(10000*.Machine$double.eps, 1-corr.u^2))*dsig2.dsig2.st

  TbT5 <- sigma2*drh.drh.st2.u
  TbT6 <- dsig2.dsig2.st*drh.drh.st.u
  TbT7 <- dsig2.dsig2.st*corr.u
  TbT8 <- (-sigma2/sqrt(pmax(10000*.Machine$double.eps, 1-corr.u^2)))*(drh.drh.st.u^2/(1-corr.u^2)+corr.u*drh.drh.st2.u)
  TbT9 <- (TbT3/sigma2)*dsig2.dsig2.st
  TbT10 <- sqrt(pmax(10000*.Machine$double.eps, 1-corr.u^2))*dsig2.dsig2.st

## Score function for (theta1, beta1, sigma1*), (theta, thet2, beta2, rho.u*, sigma2*), and rho.e*

  dl.dbe1 <- weights*q1*g1/PHI2
  dl.dbe2 <- weights*q2*g2/PHI2
  dl.drho <- weights*q1*q2*phi2/PHI2

  G1 <- c(dl.dbe1)*dat1all
  G2 <- c(dl.dbe2)*dat2all

  G1tot <- array(0,X1.d2+1)
  G2tot <- array(0,X2.d2+2)
  G3tot <- 0
  DLDP <- array(0,4) # d log-likelihood / d p, where p = (sigma1,rho_u*sigma2,sqrt(1-rho_u^2)*sigma2,rho_e)

  for (i in 1:N)
      for (j in 1:NGQ^2)
          if (uidf[i]>1) {G1tot<-G1tot+WN2[i,j]*colSums(G1[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])
                         }else{G1tot<-G1tot+WN2[i,j]*G1[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]}

  DLDP[1] <- G1tot[X1.d2+1]
  G1tot[X1.d2+1] <- G1tot[X1.d2+1]*dsig1.dsig1.st

  for (i in 1:N)
      for (j in 1:NGQ^2)
          if (uidf[i]>1) {G2tot <- G2tot+WN2[i,j]*colSums(G2[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])
                         }else{G2tot <- G2tot+WN2[i,j]*G2[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]}

  DLDP[2:3] <- G2tot[c(X2.d2+1,X2.d2+2)]
  G2tot <- c(adiag(diag(1,X2.d2),matrix(c(TbT1,TbT2,TbT3,TbT4),ncol=2))%*%matrix(G2tot))

  for (i in 1:N)
      for (j in 1:NGQ^2)
          G3tot <- G3tot+WN2[i,j]*sum(dl.drho[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])

  DLDP[4] <- G3tot
  G <- -c(G1tot,G2tot,G3tot*drh.drh.st.e)

## Observed Information Matrix
## First the negatives of second derivatives

  d2l.be1.be1<- -weights*(-w1*g1-corr.sq*phi2-g1^2/PHI2)/PHI2
  d2l.be2.be2<- -weights*(-w2*g2-corr.sq*phi2-g2^2/PHI2)/PHI2
  d2l.be1.be2<- -weights*(q1*q2)/PHI2*(phi2-g1*g2/PHI2)
  d2l.be1.rho<- -weights*(corr.sq*delta*v1-w1-g1/PHI2)*q2*phi2/PHI2
  d2l.be2.rho<- -weights*(corr.sq*delta*v2-w2-g2/PHI2)*q1*phi2/PHI2
  d2l.rho.rho<- -weights*(phi2/PHI2)*(delta^2*corr.sq*(1-delta^2*(w1^2+w2^2-2*corr.sq*w1*w2))+delta^2*w1*w2-phi2/PHI2)

  H11A <- H11B <- H11C <- matrix(0,nrow=X1.d2+1,ncol=X1.d2+1)
  H22A <- H22B <- H22C <- matrix(0,nrow=X2.d2+2,ncol=X2.d2+2)
  H33A <- H33B <- H33C <- 0
  H12A <- H12B <- H12C <- matrix(0,nrow=X1.d2+1,ncol=X2.d2+2)
  H13A <- H13B <- H13C <- matrix(0,nrow=X1.d2+1,ncol=1)
  H23A <- H23B <- H23C <- matrix(0,nrow=X2.d2+2,ncol=1)

  H1 <- c(d2l.be1.be1)*dat1all
  H2 <- c(d2l.be2.be2)*dat2all
  H12 <- c(d2l.be1.be2)*dat1all
  H13 <- c(d2l.be1.rho)*dat1all
  H23 <- c(d2l.be2.rho)*dat2all

  for (i in 1:N){
      for (j in 1:NGQ^2){
          if (uidf[i]>1) {H11A<-H11A+WN2[i,j]*crossprod(H1[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),],
          dat1all[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])}else{H11A<-H11A+WN2[i,j]*
          matrix(H1[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])%*%t(matrix(dat1all[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]))}

          if (uidf[i]>1) {H22A<-H22A+WN2[i,j]*crossprod(H2[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),],
          dat2all[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])}else{H22A<-H22A+WN2[i,j]*
          matrix(H2[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])%*%t(matrix(dat2all[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]))}

          H33A<-H33A+WN2[i,j]*sum(d2l.rho.rho[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])

          if (uidf[i]>1) {H12A<-H12A+WN2[i,j]*crossprod(H12[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),],
          dat2all[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])}else{H12A<-H12A+WN2[i,j]*
          matrix(H12[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])%*%t(matrix(dat2all[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]))}

          if (uidf[i]>1) {H13A<-H13A+WN2[i,j]*colSums(H13[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])}else{
              H13A<-H13A+WN2[i,j]*H13[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]}

          if (uidf[i]>1) {H23A<-H23A+WN2[i,j]*colSums(H23[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])}else{
              H23A<-H23A+WN2[i,j]*H23[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]}
      }
  }

  for (i in 1:N){
      Temp1<-Temp2<-Temp3<-0
      for (j in 1:NGQ^2){
          if (uidf[i]>1) {temp1<-colSums(G1[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])}else{temp1<-G1[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]}
          if (uidf[i]>1) {temp2<-colSums(G2[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])}else{temp2<-G2[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),]}
          temp3<-sum(dl.drho[seq(j+cuid[i]*NGQ^2,cuid[i+1]*NGQ^2,NGQ^2),])
          H11B<-H11B+WN2[i,j]*matrix(temp1)%*%temp1
          H22B<-H22B+WN2[i,j]*matrix(temp2)%*%temp2
          H33B<-H33B+WN2[i,j]*temp3*temp3
          H12B<-H12B+WN2[i,j]*matrix(temp1)%*%temp2
          H13B<-H13B+WN2[i,j]*matrix(temp1)*temp3
          H23B<-H23B+WN2[i,j]*matrix(temp2)*temp3
          Temp1<-Temp1+WN2[i,j]*temp1
          Temp2<-Temp2+WN2[i,j]*temp2
          Temp3<-Temp3+WN2[i,j]*temp3
      }
      H11C<-H11C+matrix(Temp1)%*%Temp1
      H22C<-H22C+matrix(Temp2)%*%Temp2
      H33C<-H33C+Temp3*Temp3
      H12C<-H12C+matrix(Temp1)%*%Temp2
      H13C<-H13C+matrix(Temp1)*Temp3
      H23C<-H23C+matrix(Temp2)*Temp3
  }

  Multi <- adiag(diag(1,X1.d2),dsig1.dsig1.st,diag(1,X2.d2),matrix(c(TbT1,TbT3,TbT2,TbT4),ncol=2),drh.drh.st.e)

  H <- t(Multi)%*% rbind(cbind(H11A-H11B+H11C, H12A-H12B+H12C, H13A-H13B+H13C),
                         cbind(t(H12A-H12B+H12C), H22A-H22B+H22C, H23A-H23B+H23C),
                         cbind(t(H13A-H13B+H13C),t(H23A-H23B+H23C),H33A-H33B+H33C)) %*%Multi

  ## Add second part of generalized chain rule

  H[X1.d2+1,X1.d2+1] <- H[X1.d2+1,X1.d2+1]-DLDP[1]*dsig1.dsig1.st
  H[X1.d2+1+X2.d2+3,X1.d2+1+X2.d2+3] <- H[X1.d2+1+X2.d2+3,X1.d2+1+X2.d2+3]-DLDP[4]*drh.drh.st2.e
  H[c(1:2)+X1.d2+1+X2.d2,c(1:2)+X1.d2+1+X2.d2] <- H[c(1:2)+X1.d2+1+X2.d2,c(1:2)+X1.d2+1+X2.d2]-DLDP[2]*matrix(c(TbT5,TbT6,TbT6,TbT7),2,2)
  H[c(1:2)+X1.d2+1+X2.d2,c(1:2)+X1.d2+1+X2.d2] <- H[c(1:2)+X1.d2+1+X2.d2,c(1:2)+X1.d2+1+X2.d2]-DLDP[3]*matrix(c(TbT8,TbT9,TbT9,TbT10),2,2)

  # So H = - d^2 logL

  ### Function returns.

  res <- -l.par
  
  
  if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0

     else{

    S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
    S <- do.call(adiag, lapply(S, unlist))

    if(l.sp1!=0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1),
                                          S[1:(X1.d2-gp1),1:(X1.d2-gp1)],
                                          0,
                      			          matrix(0,gp2,gp2),
                      			          S[(X1.d2-(gp1-1)):dim(S)[2],(X1.d2-(gp1-1)):dim(S)[2]],
                      			          0,0,0)

    if(l.sp1==0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1), 0, matrix(0,gp2,gp2), S, 0, 0, 0)
    if(l.sp1!=0 && l.sp2==0) S.h <- adiag(matrix(0,gp1,gp1), S, 0, matrix(0,gp2,gp2), 0, 0, 0)


   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
         }

         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
         H   <- H + S.h  


  p11 <- pmax( pnorm2( eta1, eta2, cov12=corr.e), 1000*.Machine$double.eps )
  p10 <- pmax( pnorm(eta1) - p11, 1000*.Machine$double.eps )
  p01 <- pmax( pnorm(eta2) - p11, 1000*.Machine$double.eps )
  p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )


         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, 
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho) 
  
}

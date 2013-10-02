bprobgHs.NRE <- function(params, BivD, nC, nu,  gevOut, xi, H.n=NULL, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1=NULL, X2=NULL, weights=weights, X1.d2, X2.d2, pPen1=NULL, pPen2=NULL, sp=NULL, qu.mag=NULL, gp1, gp2, fp, l.sp1, l.sp2, K=NULL, n,      N,      cuid,      uidf,      masses=NULL, NGQ,      dat1all, dat2all,           W){

## Papameters (transformed or not), linear predictors, and some other definitions

    sigma1.st <- params[X1.d2+1]
    corr.st.u <- params[X1.d2+X2.d2+2]
    sigma2.st <- params[X1.d2+X2.d2+3]
    corr.st.e <- params[X1.d2+X2.d2+4]

    sigma1 <- exp(sigma1.st)
    corr.u <- tanh(corr.st.u)
    sigma2 <- exp(sigma2.st)

    epsilon <- .Machine$double.eps*10^6
    if(BivD %in% c("N","T")      ){corr.e <- tanh(corr.st.e); if(corr.e %in% c(-1,1)) corr.e <- sign(corr.e)*0.9999999}
    if(BivD=="F")                  corr.e <- corr.st.e + epsilon
    if(BivD %in% c("C0", "C180") ) corr.e <- exp(corr.st.e) + epsilon
    if(BivD %in% c("C90","C270") ) corr.e <- -( exp(corr.st.e) + epsilon )
    if(BivD %in% c("J0", "J180") ) corr.e <- exp(corr.st.e) + 1 + epsilon
    if(BivD %in% c("J90","J270") ) corr.e <- -( exp(corr.st.e) + 1 + epsilon )
    if(BivD %in% c("G0", "G180") ) corr.e <- exp(corr.st.e) + 1
    if(BivD %in% c("G90","G270") ) corr.e <- -( exp(corr.st.e) + 1 )

    eta1 <- dat1all%*%c(params[1:X1.d2],sigma1)
    eta2 <- dat2all%*%c(params[(X1.d2+2):(X1.d2+X2.d2+1)],corr.u*sigma2,sqrt(pmax(10000*.Machine$double.eps, 1-corr.u^2))*sigma2)

    p1 <- pnorm(eta1)
    p2 <- pnorm(eta2)
    d.n1   <- dnorm(eta1)
    d.n2   <- dnorm(eta2)

## Leading to log-likelihood (l.par) and normalized weights (WN2)

    if(BivD=="N") C.copula <- pmax( abs(pnorm2( eta1, eta2, cov12=corr.e)), 1000*.Machine$double.eps ) else C.copula <- BiCopCDF(p1,p2, nC, par=corr.e, par2=nu)

    p11 <- pmax( C.copula, 1000*.Machine$double.eps )
    p10 <- pmax( p1 - p11, 1000*.Machine$double.eps )
    p01 <- pmax( p2 - p11, 1000*.Machine$double.eps )
    p00 <- pmax( 1- p11 - p10 - p01, 1000*.Machine$double.eps )

    l.par1 <- weights*( y1.y2*log(p11)+y1.cy2*log(p10)+cy1.y2*log(p01)+cy1.cy2*log(p00) )

    WN <- matrix(0,nrow=N,ncol=NGQ^2)
    for (r in 1:N)
        for (c in 1:NGQ^2)
            WN[r,c] <- exp(sum(l.par1[seq(c+cuid[r]*NGQ^2,cuid[r+1]*NGQ^2,NGQ^2)]))*W[c]

    l.par2 <- log(apply(WN,1,sum)) 
    l.par  <- sum(l.par2)

    WN2 <- (1/apply(WN,1,sum))*WN

## Derivatives for chain rule

    drh.drh.st.e <- 4*exp(2*corr.st.e)/(exp(2*corr.st.e)+1)^2 # <- (1/cosh(corr.st.e)^2)
    drh.drh.st.u <- 4*exp(2*corr.st.u)/(exp(2*corr.st.u)+1)^2 # <- (1/cosh(corr.st.u)^2)
    dsig1.dsig1.st <- exp(sigma1.st)
    dsig2.dsig2.st <- exp(sigma2.st)
    drh.drh.st2.e <- 8*exp(2*corr.st.e)*(1-exp(2*corr.st.e))/(exp(2*corr.st.e)+1)^3
    drh.drh.st2.u <- 8*exp(2*corr.st.u)*(1-exp(2*corr.st.u))/(exp(2*corr.st.u)+1)^3

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

## Derivatives from Ros' files

    dH <- copgHs(p1,p2,corr.e,corr.st.e,BivD,nC,nu)

    c.copula.be1   <- dH$c.copula.be1
    c.copula.be2   <- dH$c.copula.be2
    c.copula.theta <- dH$c.copula.theta

    c.copula2.be1 <- dH$c.copula2.be1
    bit1.b1b1 <- c.copula2.be1*(d.n1)^2-c.copula.be1*d.n1*eta1
    bit2.b1b1 <- -d.n1*eta1-bit1.b1b1
    bit3.b1b1 <- -bit1.b1b1
    bit4.b1b1 <- -bit2.b1b1

    c.copula2.be2 <- dH$c.copula2.be2
    bit1.b2b2 <- c.copula2.be2*(d.n2)^2-c.copula.be2*d.n2*eta2
    bit2.b2b2 <- -bit1.b2b2
    bit3.b2b2 <- -d.n2*eta2-bit1.b2b2
    bit4.b2b2 <- -bit3.b2b2

    c.copula2.be1be2 <- dH$c.copula2.be1be2
    bit1.b1b2 <- c.copula2.be1be2 * d.n1 *d.n2
    bit2.b1b2 <- -bit1.b1b2
    bit3.b1b2 <- -bit1.b1b2
    bit4.b1b2 <- bit1.b1b2

    c.copula2.be1th <- dH$c.copula2.be1th
    bit1.b1th <- c.copula2.be1th*d.n1
    bit2.b1th <- -bit1.b1th
    bit3.b1th <- -bit1.b1th
    bit4.b1th <- bit1.b1th

    c.copula2.be2th <- dH$c.copula2.be2th
    bit1.b2th <- c.copula2.be2th*d.n2
    bit2.b2th <- -bit1.b2th
    bit3.b2th <- -bit1.b2th
    bit4.b2th <- bit1.b2th

    bit1.th2 <- dH$bit1.th2
    bit2.th2 <- -bit1.th2
    bit3.th2 <- -bit1.th2
    bit4.th2 <- bit1.th2

    dl.dbe1 <-  weights*d.n1*( (y1.y2*c.copula.be1/p11)  +
                      (y1.cy2*(1-c.copula.be1)/p10) +
                      (cy1.y2*c.copula.be1/(-p01)) +
                      (cy1.cy2*(c.copula.be1-1)/p00) )

    dl.dbe2 <-  weights*d.n2*( (y1.y2*c.copula.be2/p11)  +
                            (y1.cy2*c.copula.be2/(-p10)) +
                                (cy1.y2*(1-c.copula.be2)/(p01)) +
                                (cy1.cy2*(c.copula.be2-1)/p00) )

    dl.drho <- weights*( y1.y2*c.copula.theta/p11+y1.cy2*(-c.copula.theta)/p10 +
                       cy1.y2*(-c.copula.theta)/p01+cy1.cy2*c.copula.theta/p00 )

    d2l.be1.be1  <- -weights*(y1.y2*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2+
                              y1.cy2*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2+
                              cy1.y2*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2+
                              cy1.cy2*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2 )

    d2l.be2.be2  <- -weights*(y1.y2*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              y1.cy2*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2+
                              cy1.y2*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2+
                              cy1.cy2*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2 )

    d2l.be1.be2  <- -weights*(y1.y2*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2+
                              y1.cy2*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2+
                              cy1.y2*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2+
                              cy1.cy2*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2 )

    d2l.be1.rho  <- -weights*(y1.y2*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )

    d2l.be2.rho  <- -weights*(y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                              cy1.y2*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              cy1.cy2*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )

    d2l.rho.rho  <- -weights*(y1.y2*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                              y1.cy2*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2+
                              cy1.y2*(bit3.th2*p01-(-c.copula.theta)^2)/p01^2+
                              cy1.cy2*(bit4.th2*p00-c.copula.theta^2)/p00^2 )

## Score function for (theta1, beta1, sigma1*), (theta, thet2, beta2, rho.u*, sigma2*), and rho.e*

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
    G <- -c(G1tot,G2tot,G3tot)

## Observed Information Matrix

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

  Multi <- adiag(diag(1,X1.d2),dsig1.dsig1.st,diag(1,X2.d2),matrix(c(TbT1,TbT3,TbT2,TbT4),ncol=2),1)

  H <- t(Multi)%*% rbind(cbind(H11A-H11B+H11C, H12A-H12B+H12C, H13A-H13B+H13C),
                         cbind(t(H12A-H12B+H12C), H22A-H22B+H22C, H23A-H23B+H23C),
                         cbind(t(H13A-H13B+H13C),t(H23A-H23B+H23C),H33A-H33B+H33C)) %*%Multi

  ## Add second part of generalized chain rule

  H[X1.d2+1,X1.d2+1] <- H[X1.d2+1,X1.d2+1]-DLDP[1]*dsig1.dsig1.st
  H[X1.d2+1+X2.d2+3,X1.d2+1+X2.d2+3] <- H[X1.d2+1+X2.d2+3,X1.d2+1+X2.d2+3]-DLDP[4]*0
  H[c(1:2)+X1.d2+1+X2.d2,c(1:2)+X1.d2+1+X2.d2] <- H[c(1:2)+X1.d2+1+X2.d2,c(1:2)+X1.d2+1+X2.d2]-DLDP[2]*matrix(c(TbT5,TbT6,TbT6,TbT7),2,2)
  H[c(1:2)+X1.d2+1+X2.d2,c(1:2)+X1.d2+1+X2.d2] <- H[c(1:2)+X1.d2+1+X2.d2,c(1:2)+X1.d2+1+X2.d2]-DLDP[3]*matrix(c(TbT8,TbT9,TbT9,TbT10),2,2)

  # So H = - d^2 logL

  ### Function returns.

  res <- -l.par
     
if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0

     else{
        
    dimP1 <- dimP2 <- 0     
    S1 <- S2 <- matrix(0,1,1)   

    S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
    S <- do.call(adiag, lapply(S, unlist))

    ma1 <- matrix(0,gp1,gp1) 
    if(length(pPen1)!=0){ indP1 <- qu.mag$off[1]:(qu.mag$off[1]+qu.mag$rank[1]-1)
                          dimP1 <- length(indP1)
                          ma1[indP1,indP1] <- S[1:dimP1,1:dimP1]
                                } 
    ma2 <- matrix(0,gp2,gp2)
    if(length(pPen2)!=0){ 
                          indP2 <- (qu.mag$off[l.sp1+1]-X1.d2):(-X1.d2+qu.mag$off[l.sp1+1]+qu.mag$rank[l.sp1+1]-1)
                          dimP2 <- length(indP2)
                          ma2[indP2,indP2] <- S[indP2+X1.d2,indP2+X1.d2]
                                }                                 
    
    lP1 <- length(pPen1); lP2 <- length(pPen2) 
    
    if((lP1!=0 && l.sp1>1) || (lP1==0 && l.sp1>0)) S1 <- S[(dimP1+1):(dimP1+X1.d2-gp1),(dimP1+1):(dimP1+X1.d2-gp1)]
    if((lP2!=0 && l.sp2>1) || (lP2==0 && l.sp2>0)){dS1 <- dim(S1)[2]; if(dS1==1) dS1 <- 0; 
                                                   S2 <- S[(dimP1+dimP2+dS1+1):dim(S)[2],(dimP1+dimP2+dS1+1):dim(S)[2]]}
    
    lS1 <- length(S1); lS2 <- length(S2) 
    
    if(lS1==1 && lS2==1) S.h <- adiag(ma1, 0, ma2, 0, 0, 0)
    if(lS1 >1 && lS2==1) S.h <- adiag(ma1, S1, 0, ma2, 0, 0, 0)
    if(lS1==1 && lS2 >1) S.h <- adiag(ma1, 0, ma2, S2, 0 , 0, 0)
    if(lS1 >1 && lS2 >1) S.h <- adiag(ma1, S1, 0, ma2, S2, 0, 0, 0)
        
   
   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
   
         }         
         
         
         S.res <- res
         res <- S.res + S.h1
         G   <- G + S.h2
         H   <- H + S.h

         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, # l.par=l.par2, # do I need here l.par2 or l.par1?
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2,
              d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho)
}

EBNormal<-function(params,y1,y2,X1,X2,X1.d2,X2.d2,N,cuid,uidf,tol,M){

    ## Papameters and responses

    sigma1.st <- params[X1.d2+1]
    corr.st.u <- params[X1.d2+X2.d2+2]
    sigma2.st <- params[X1.d2+X2.d2+3]
    corr.st.e <- params[X1.d2+X2.d2+4]
    
    sigma1 <- exp(sigma1.st)
    corr.u <- tanh(corr.st.u)
    sigma2 <- exp(sigma2.st)
    corr.e <- tanh(corr.st.e)
    
    Sigma <- matrix(c(sigma1^2,sigma1*sigma2*corr.u,sigma1*sigma2*corr.u,sigma2^2),nrow=2,ncol=2)
    
    q1ALL <- 2*y1-1
    q2ALL <- 2*y2-1
    
    eta1ALL <- X1%*%matrix(params[1:X1.d2]) # want to exclude RE parameters, why?
    eta2ALL <- X2%*%matrix(params[(X1.d2+2):(X1.d2+X2.d2+1)]) # want to exclude RE parameters, why?

    ## Calculate and store maximized likelihood wrt random effects

    supF <- array(0,N)

    ## Subject specific analysis

    for (i in 1:N){
        #if ((i/20)==round(i/20)) {cat('..')}
        u <- matrix(0,nrow=2)

        q1 <- q1ALL[(cuid[i]+1):(cuid[i+1])]
        q2 <- q2ALL[(cuid[i]+1):(cuid[i+1])]
        corr.sq <- q1*q2*corr.e
        delta <- 1/sqrt(pmax(10000*.Machine$double.eps,1-corr.sq^2))

        G <- matrix(10,nrow=2)

        while(max(abs(G)) > tol && sum(as.numeric(G=="NaN"))==0){

            eta1 <- eta1ALL[(cuid[i]+1):(cuid[i+1])]+u[1]
            eta2 <- eta2ALL[(cuid[i]+1):(cuid[i+1])]+u[2]
            w1 <- q1*eta1
            w2 <- q2*eta2
            v1 <- delta*(w2-corr.sq*w1)
            v2 <- delta*(w1-corr.sq*w2)
            g1 <- dnorm(w1)*pnorm(v1)
            g2 <- dnorm(w2)*pnorm(v2)
            phi2 <- dnorm2(w1,w2,rho=corr.sq)
            PHI2 <- pmax(abs(pnorm2(w1,w2,cov12=corr.sq)),1000*.Machine$double.eps)

            ## Score function for ui1 and ui2

            dl.dbe1 <- q1*g1/PHI2
            dl.dbe2 <- q2*g2/PHI2
            G <- matrix(c(sum(dl.dbe1),sum(dl.dbe2)))

            ## Observed Information matrix for the (random) intercepts: H = - d^2 logL

            d2l.be1.be1<- -(-w1*g1-corr.sq*phi2-g1^2/PHI2)/PHI2
            d2l.be2.be2<- -(-w2*g2-corr.sq*phi2-g2^2/PHI2)/PHI2
            d2l.be1.be2<- -(q1*q2)/PHI2*(phi2-g1*g2/PHI2)
            H <- matrix(c(sum(d2l.be1.be1),sum(d2l.be1.be2),sum(d2l.be1.be2),sum(d2l.be2.be2)),ncol=2)

            u <- u + ginv(H)%*%G
        }
        eta1 <- eta1ALL[(cuid[i]+1):(cuid[i+1])]+u[1]
        eta2 <- eta2ALL[(cuid[i]+1):(cuid[i+1])]+u[2]
        w1 <- q1*eta1
        w2 <- q2*eta2
        PHI2 <- pmax(abs(pnorm2(w1,w2,cov12=corr.sq)),1000*.Machine$double.eps)
        supF[i] <- prod(PHI2)
    }

    ## Posterior samples

    EB1 <- matrix(0,nrow=N,ncol=M)
    EB2 <- matrix(0,nrow=N,ncol=M)

    for (i in 1:N){
        #if ((i/10)==round(i/10)) {cat('..')}
        sample <- 0
        q1 <- q1ALL[(cuid[i]+1):(cuid[i+1])]
        q2 <- q2ALL[(cuid[i]+1):(cuid[i+1])]
        corr.sq <- q1*q2*corr.e
        while(sample<M){
            up<-rmvnorm(1,c(0,0),Sigma)
            eta1<-eta1ALL[(cuid[i]+1):(cuid[i+1])]+up[1]
            eta2<-eta2ALL[(cuid[i]+1):(cuid[i+1])]+up[2]
            w1<-q1*eta1
            w2<-q2*eta2
            PHI2<-pmax(abs(pnorm2(w1,w2,cov12=corr.sq)),1000*.Machine$double.eps)
            if (runif(1) < prod(PHI2)/supF[i]){
                sample <- sample+1
                EB1[i,sample] <- up[1]
                EB2[i,sample] <- up[2]
            }
        }
    }
    
    eb.u1 <- apply(EB1, 1, mean) # should this return the same sd as the estimated one? 
    eb.u2 <- apply(EB2, 1, mean)
    Eb.u1 <- rep(eb.u1,uidf)
    Eb.u2 <- rep(eb.u2,uidf)

    ## Function returns:

    list(Eb.u1=Eb.u1,Eb.u2=Eb.u2,eb.u1=eb.u1,eb.u2=eb.u2)
    
}

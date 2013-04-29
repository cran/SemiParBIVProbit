NP.qu <- function(SemiParFit, y1, y2, q1, q2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, X1.d2, X2.d2, qu.mag=NULL, gp1, gp2, 
                    fp, l.sp1, l.sp2, weights, K, n, N, cuid, uidf, NGQ=NULL, dat1all=NULL, dat2all=NULL, W=NULL){
                            
        T.sv <- bprobNP.H(SemiParFit$fit$argument, y1=y1, y2=y2, 
                         q1=q1, q2=q2, y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1, 
                         X1=X1, X2=X2, 
                         X1.d2=X1.d2, X2.d2=X2.d2, sp=SemiParFit$sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, weights=weights,
                         K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=SemiParFit$masses, NGQ=NGQ, dat1all=dat1all, dat2all=dat2all, W=W)

        npar <- K + X1.d2 + K + X2.d2 + 1 + K - 1
        H.cor <- G.cor <- matrix(0,nrow=npar,ncol=npar)

        for (i in 1:N){
              h <- matrix(0,nrow=npar,ncol=1)
           for (l in 1:K){

              c1 <- T.sv$dl.dbe1[(cuid[i]+1):(cuid[i+1]),l]*X1(l)[(cuid[i]+1):(cuid[i+1]),]
              if (uidf[i]>1) c1<-apply(c1,2,sum)

              c2 <- T.sv$dl.dbe2[(cuid[i]+1):(cuid[i+1]),l]*X2(l)[(cuid[i]+1):(cuid[i+1]),]
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

        u1 <- SemiParFit$fit$argument[1:K]
        u2 <- SemiParFit$fit$argument[(K+X1.d2+1):(K+X1.d2+K)]

        nw <- T.sv$Wp3/matrix(rep(apply(T.sv$Wp3,1,sum),K),ncol=K)

        eb.u1 <- apply(t(u1*t(nw)),1,sum)
        eb.u2 <- apply(t(u2*t(nw)),1,sum)

        Eb.u1 <- rep(eb.u1,uidf)
        Eb.u2 <- rep(eb.u2,uidf)

    list(He=He,logL=logL,eb.u1=eb.u1,eb.u2=eb.u2,Eb.u1=Eb.u1,Eb.u2=Eb.u2,T.sv=T.sv)

    }


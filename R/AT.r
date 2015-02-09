AT <- function(x, eq, nm.bin = "", E = TRUE, treat = TRUE, delta = FALSE, prob.lev = 0.05, s.meth = "svd", n.sim = 500, naive = FALSE){

stopm <- "Calculation of this average treatment effect is valid for recursive binary models only."
end <- 0

if(x$ig[[1]]$response %in% x$ig[[2]]$pred.names ) end <- 1
if(x$ig[[2]]$response %in% x$ig[[1]]$pred.names ) end <- 1
if(x$Model=="BSS" || x$Model=="BPO" || end==0) stop(stopm)

if(is.character(nm.bin)==FALSE) stop("nm.bin is not a character!")

if(x$PL!="P") delta <- FALSE

est.ATb <- NA
ind <- list()

good <- x$fit$good
epsilon <- sqrt(.Machine$double.eps)







if(naive==FALSE){

p.rho <- length(coef(x))

if(x$PL!="P" && x$fitPL!="fixed"){if(x$eqPL=="both") p.rho <- p.rho - 2 else p.rho <- p.rho - 1} 

ind[[1]] <- 1:x$X1.d2 
ind[[2]] <- x$X1.d2+(1:x$X2.d2)

}




if(eq==1){ X.int <- x$X1[good,]
    if(naive==FALSE){
           X.noi <- x$X2[good,]
           ind.int <- ind[[1]]
           ind.noi <- ind[[2]] 
           etap.noi <- x$eta2 }
                 
}

if(eq==2){ X.int <- x$X2[good,]
    if(naive==FALSE){
           X.noi <- x$X1[good,]
           ind.int <- ind[[2]]
           ind.noi <- ind[[1]] 
           etap.noi  <- x$eta1}
            
}


if(naive==FALSE){              
coef.int <- as.numeric(coef(x)[ind.int])
coef.noi <- as.numeric(coef(x)[ind.noi])        
}





X.int <- X.int[good,]       
if(naive==FALSE) X.noi <- X.noi[good,]   
       

                      
d0 <- d1 <- X.int
d0[,nm.bin] <- 0
d1[,nm.bin] <- 1



if(naive==FALSE){
eti1 <- d1%*%coef.int 
eti0 <- d0%*%coef.int 
etno <- etap.noi      
}



if(E==TRUE) ind.excl <- rep(TRUE,x$n) else{ 

      if(treat==TRUE) ind.excl <- as.logical(X.int[,nm.bin]) 
                 else ind.excl <- as.logical(X.int[,nm.bin])!=TRUE
      
    } 






###############
# Calculate ATE
###############

if(naive==FALSE){





if(x$PL=="P"){

p.int1 <- pnorm(eti1)
p.int0 <- pnorm(eti0)


d.int1 <- dnorm(eti1)
d.int0 <- dnorm(eti0)

p.etn  <- pnorm(etno)
d.etn  <- dnorm(etno)
}




if(x$PL=="PP"){

if(eq==1){

p.int1 <- pnorm(eti1)^x$xi1
p.int0 <- pnorm(eti0)^x$xi1

d.int1 <- pnorm(eti1)^(x$xi1 - 1) * (x$xi1 * dnorm(eti1))
d.int0 <- pnorm(eti0)^(x$xi1 - 1) * (x$xi1 * dnorm(eti0))

p.etn  <- pnorm(etno)^x$xi2   
d.etn  <- pnorm(etno)^(x$xi2 - 1) * (x$xi2 * dnorm(etno))

         }   

if(eq==2){

p.int1 <- pnorm(eti1)^x$xi2
p.int0 <- pnorm(eti0)^x$xi2

d.int1 <- pnorm(eti1)^(x$xi2 - 1) * (x$xi2 * dnorm(eti1))
d.int0 <- pnorm(eti0)^(x$xi2 - 1) * (x$xi2 * dnorm(eti0))

p.etn  <- pnorm(etno)^x$xi1   
d.etn  <- pnorm(etno)^(x$xi1 - 1) * (x$xi1 * dnorm(etno))
        }

}



if(x$PL=="RPP"){

if(eq==1){

p.int1 <- 1-pnorm(-eti1)^x$xi1
p.int0 <- 1-pnorm(-eti0)^x$xi1

d.int1 <- pnorm(-eti1)^(x$xi1 - 1) * (x$xi1 * dnorm(-eti1))
d.int0 <- pnorm(-eti0)^(x$xi1 - 1) * (x$xi1 * dnorm(-eti0))

p.etn  <- 1-pnorm(-etno)^x$xi2   
d.etn  <- pnorm(-etno)^(x$xi2 - 1) * (x$xi2 * dnorm(-etno))

         }   

if(eq==2){

p.int1 <- 1-pnorm(-eti1)^x$xi2
p.int0 <- 1-pnorm(-eti0)^x$xi2

d.int1 <- pnorm(-eti1)^(x$xi2 - 1) * (x$xi2 * dnorm(-eti1))
d.int0 <- pnorm(-eti0)^(x$xi2 - 1) * (x$xi2 * dnorm(-eti0))

p.etn  <- 1-pnorm(-etno)^x$xi1   
d.etn  <- pnorm(-etno)^(x$xi1 - 1) * (x$xi1 * dnorm(-etno))
        }

}


  if(x$PL=="SN"){  
  
  del1 <- -x$xi1/sqrt(1+x$xi1^2)
  del2 <- -x$xi2/sqrt(1+x$xi2^2)


  if(eq==1){

p.int1 <- 2*pbinorm( eti1, 0, cov12=del1)
p.int0 <- 2*pbinorm( eti0, 0, cov12=del1)

d.int1 <- 2*dnorm(eti1)*pnorm(x$xi1*eti1)
d.int0 <- 2*dnorm(eti0)*pnorm(x$xi1*eti0)

p.etn  <- 2*pbinorm( etno, 0, cov12=del2)
d.etn  <- 2*dnorm(etno)*pnorm(x$xi2*etno)

         }   

if(eq==2){

p.int1 <- 2*pbinorm( eti1, 0, cov12=del2)
p.int0 <- 2*pbinorm( eti0, 0, cov12=del2)

d.int1 <- 2*dnorm(eti1)*pnorm(x$xi2*eti1)
d.int0 <- 2*dnorm(eti0)*pnorm(x$xi2*eti0)

p.etn  <- 2*pbinorm( etno, 0, cov12=del1)
d.etn  <- 2*dnorm(etno)*pnorm(x$xi1*etno)
        }    
      
  
   
  }
  




if( is.null(x$X3) ) {
  if(x$BivD == "N" ) {ass.p <- x$rho; ass.pst <- coef(x)["athrho"] } else{ ass.p <- x$theta; ass.pst <- coef(x)["theta.star"]}  
}

if( !is.null(x$X3) ) {
  if(x$BivD == "N" ) {ass.p <- x$rho; ass.pst <- x$fit$theta.star } else{ ass.p <- x$theta; ass.pst <- x$fit$theta.star}  
}


p.int1 <- pmax(p.int1, epsilon ) 
p.int1 <- ifelse(p.int1==1, 0.9999999, p.int1)
p.int0 <- pmax(p.int0, epsilon )
p.int0 <- ifelse(p.int0==1, 0.9999999,p.int0)


p.etn  <- pmax(p.etn, epsilon )
p.etn <- ifelse(p.etn==1, 0.9999999,p.etn)


   C.11  <- pmax( BiCDF(p.int1, p.etn, x$nC, ass.p) , epsilon )
   C.10  <- p.int0 - pmax( BiCDF(p.int0, p.etn, x$nC, ass.p) , epsilon )
                                 

   est.AT <- mean(   (C.11/p.etn - C.10/(1-p.etn))[ind.excl],na.rm = TRUE   )

} # end naive condition 




if(naive == TRUE){

if(eq==1) ngam <- x$gam1
if(eq==2) ngam <- x$gam2

eti1 <- d1%*%coef(ngam) 
eti0 <- d0%*%coef(ngam) 

p.int1 <- pmax(pnorm(eti1), epsilon ) 
p.int1 <- ifelse(p.int1==1,0.9999999,p.int1) 
p.int0 <- pmax(pnorm(eti0), epsilon ) 
p.int0 <- ifelse(p.int0==1,0.9999999,p.int0) 

est.AT <- mean(  (p.int1 - p.int0)[ind.excl] , na.rm = TRUE    ) 


}


###############





      
if(delta==TRUE){

  if(naive == FALSE){

    if(x$BivD %in% c("N")      )     add.b <- 1/cosh(ass.pst)^2
    if(x$BivD=="F")                  add.b <- 1
    if(x$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) add.b <-  exp(ass.pst)
    if(x$BivD %in% c("C90","C270","J90","J270","G90","G270") ) add.b <- -exp(ass.pst)
   
   dC1 <- copgHs(p1=p.etn,p2=p.int1,teta=ass.p,teta.st=ass.pst,BivD=x$BivD,nC=x$nC,nu=NULL,xi1=NULL,xi2=NULL,eta1=etno,eta2=eti1,PL=x$PL,eqPL=x$eqPL)
   dC0 <- copgHs(p1=p.etn,p2=p.int0,teta=ass.p,teta.st=ass.pst,BivD=x$BivD,nC=x$nC,nu=NULL,xi1=NULL,xi2=NULL,eta1=etno,eta2=eti0,PL=x$PL,eqPL=x$eqPL)

   dATT.noint <- ( (dC1$c.copula.be1*p.etn-C.11)/p.etn^2 + (dC0$c.copula.be1*(1-p.etn)-C.10)/(1-p.etn)^2)*d.etn 
   dATT.noint <- colMeans( (c(dATT.noint)*X.noi)[ind.excl,] ) 

   dATT.int   <- colMeans( (c( (dC1$c.copula.be2*d.int1)/p.etn )*d1)[ind.excl,] ) - colMeans( (c( (1-dC0$c.copula.be2)*d.int0/(1-p.etn) )*d0)[ind.excl,])

   if( is.null(x$X3) ) dATT.tet   <- mean(((dC1$c.copula.theta/add.b)/p.etn + (dC0$c.copula.theta/add.b)/(1-p.etn))[ind.excl])  
   
   if( !is.null(x$X3) ) dATT.tet <- colMeans( (c( (dC1$c.copula.theta/add.b)/p.etn + (dC0$c.copula.theta/add.b)/(1-p.etn) )*x$X3[good,])[ind.excl,]  )  


   if(eq==2) dATT <- c(dATT.noint,dATT.int,dATT.tet) else dATT <- c(dATT.int,dATT.noint,dATT.tet) 


   var <- bprobgHs(params=coef(x), PL=x$PL, eqPL=x$eqPL, respvec=x$respvec, VC=x$VC, 
                   sp=x$sp, qu.mag=x$qu.mag, AT=TRUE)$hessian


   var.eig <- eigen(var, symmetric=TRUE)                    
   if(min(var.eig$values) < epsilon) var.eig$values[which(var.eig$values < epsilon)] <- 0.0000001
   var <- var.eig$vectors%*%tcrossprod(diag(1/var.eig$values),var.eig$vectors)  

   delta.AT <- sqrt( t(dATT)%*%var%*%dATT )
   
   }
   
   
   if(naive == TRUE){
   
   var <- ngam$Vp   
   dATT <- colMeans(  ( c(dnorm(eti1))*d1 - c(dnorm(eti0))*d0 )[ind.excl,]  )
   delta.AT <- sqrt( t(dATT)%*%var%*%dATT )
       
   }
   
                      
}



if(delta==FALSE){


if(naive == TRUE){ 


bs <- rmvnorm(n.sim, mean = coef(ngam), sigma=ngam$Vp, method=s.meth) 

 peti1s <- pmax(pnorm(d1%*%t(bs)), epsilon )  
 peti1s <- ifelse(peti1s==1,0.9999999,peti1s)  
 peti0s <- pmax(pnorm(d0%*%t(bs)), epsilon ) 
 peti0s <- ifelse(peti0s==1,0.9999999,peti0s)  


 est.ATb <- colMeans(  (peti1s - peti0s)[ind.excl,] , na.rm = TRUE    )


}



if(naive == FALSE){ 

bs <- rmvnorm(n.sim, mean = coef(x), sigma=x$Vb, method=s.meth)

 eti1s <- d1%*%t(bs[,ind.int])    
 eti0s <- d0%*%t(bs[,ind.int])   
 etnos <- X.noi%*%t(bs[,ind.noi]) 
 
                          
if(x$PL=="P"){
	p.int1s <- pnorm(eti1s)
	p.int0s <- pnorm(eti0s)
	p.etns  <- pnorm(etnos)

             }


if(x$PL=="PP"){                    
 if(eq==1){
  p.int1s <- pnorm(eti1s)^x$xi1   
  p.int0s <- pnorm(eti0s)^x$xi1
  p.etns  <- pnorm(etnos)^x$xi2 
           }   
 if(eq==2){
  p.int1s <- pnorm(eti1s)^x$xi2
  p.int0s <- pnorm(eti0s)^x$xi2
  p.etns  <- pnorm(etnos)^x$xi1 
           }
}


if(x$PL=="RPP"){
 if(eq==1){
  p.int1s <- 1-pnorm(-eti1s)^x$xi1
  p.int0s <- 1-pnorm(-eti0s)^x$xi1
  p.etns  <- 1-pnorm(-etnos)^x$xi2  
          }   
 if(eq==2){
  p.int1s <- 1-pnorm(-eti1s)^x$xi2
  p.int0s <- 1-pnorm(-eti0s)^x$xi2
  p.etns  <- 1-pnorm(-etnos)^x$xi1  
          }
}


  
if(x$PL=="SN"){  

dim1 <- dim(eti1s)[1]
dim2 <- dim(eti1s)[2]
      
  if(eq==1){

p.int1s <- matrix(2*pbinorm( as.vector(eti1s), 0, cov12=del1),dim1,dim2)
p.int0s <- matrix(2*pbinorm( as.vector(eti0s), 0, cov12=del1),dim1,dim2)
p.etns  <- matrix(2*pbinorm( as.vector(etnos), 0, cov12=del2),dim1,dim2)

         }   

if(eq==2){

p.int1s <- matrix(2*pbinorm( as.vector(eti1s), 0, cov12=del2),dim1,dim2)
p.int0s <- matrix(2*pbinorm( as.vector(eti0s), 0, cov12=del2),dim1,dim2)
p.etns  <- matrix(2*pbinorm( as.vector(etnos), 0, cov12=del1),dim1,dim2)

        }    
   
  }  
  
  
  
  
  
 if( !is.null(x$X3) ) etds <- x$X3[good,]%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)])
 if(  is.null(x$X3) ) etds <- bs[, p.rho]  
  
  
#if(x$BivD!="N"){

   if(x$BivD=="F")                   ass.ps <- etds + epsilon
   if(x$BivD %in% c("N"))           {ass.ps <- tanh(etds); ass.ps <- ifelse(ass.ps == -1, -0.9999999, ass.ps) 
                                                           ass.ps <- ifelse(ass.ps ==  1,  0.9999999, ass.ps) }
  
   if(x$BivD %in% c("C0", "C180") )  ass.ps <-   exp(etds) + epsilon
   if(x$BivD %in% c("C90","C270") )  ass.ps <- -(exp(etds) + epsilon)

   if(x$BivD %in% c("J0", "J180") )  ass.ps <-    1 + exp(etds) + epsilon
   if(x$BivD %in% c("J90","J270") )  ass.ps <- -( 1 + exp(etds) + epsilon)
 
   if(x$BivD %in% c("G0", "G180") )  ass.ps <-    1 + exp(etds)
   if(x$BivD %in% c("G90","G270") )  ass.ps <- -( 1 + exp(etds) )
   
   ass.ps <- ifelse(ass.ps == Inf ,  8.218407e+307, ass.ps )
   ass.ps <- ifelse(ass.ps == -Inf, -8.218407e+307, ass.ps )

#}




p.int1s <- pmax(p.int1s, epsilon )
p.int1s <- ifelse(p.int1s==1,0.9999999,p.int1s)   
p.int0s <- pmax(p.int0s, epsilon )
p.int0s <- ifelse(p.int0s==1,0.9999999,p.int0s) 
p.etns  <- pmax(p.etns, epsilon )
p.etns  <- ifelse(p.etns==1,0.9999999,p.etns) 






if( is.null(x$X3) ) ass.ps <- t(matrix(ass.ps)) 


for(i in 1:n.sim){ 

 C.11 <- pmax( BiCDF(p.int1s[,i], p.etns[,i], x$nC, ass.ps[,i]) , epsilon )
 C.10 <- p.int0s[,i] - pmax( BiCDF(p.int0s[,i], p.etns[,i], x$nC, ass.ps[,i]) , epsilon )

 est.ATb[i] <- mean(   (C.11/p.etns[,i] - C.10/(1-p.etns[,i]))[ind.excl],na.rm = TRUE   )
                  }
                  

               
}  # end big loop   


}


           
if(delta==TRUE){esd.AT <- delta.AT*qnorm(prob.lev/2,lower.tail = FALSE) 
                   CIs <- c(est.AT - esd.AT, est.AT + esd.AT)
               }else CIs <- as.numeric(quantile(est.ATb,c(prob.lev/2,1-prob.lev/2),na.rm = TRUE))

res <- c(CIs[1],est.AT,CIs[2])
out <- list(res=res, prob.lev=prob.lev, est.ATb=est.ATb)
 
class(out) <- "AT"

out

}




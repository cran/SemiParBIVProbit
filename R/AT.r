AT <- function(x, eq, nm.bin = "", E = TRUE, treat = TRUE, delta = FALSE, prob.lev = 0.05, s.meth = "svd", n.sim = 1000){

stopm <- "Calculation of this average treatment effect is valid for recursive binary models only."
end <- 0
if(x$ig[[1]]$response %in% x$ig[[2]]$pred.names ) end <- 1
if(x$ig[[2]]$response %in% x$ig[[1]]$pred.names ) end <- 1
if(x$Model=="BSS" || x$Model=="BPO" || end==0) stop(stopm)

if(is.character(nm.bin)==FALSE) stop("nm.bin is not a character!")

# delta method does not allow the use of the fisher information

if(x$PL!="P") delta <- FALSE

est.ATb <- NA
ass.ps2 <- ass.p2 <- NULL
ind <- list()

p.rho <- length(coef(x))

if(x$PL!="P" && x$fitPL!="fixed"){if(x$eqPL=="both") p.rho <- p.rho - 2 else p.rho <- p.rho - 1} 

good <- x$fit$good
epsilon <- .Machine$double.eps*10^6

ind[[1]] <- 1:x$X1.d2 
ind[[2]] <- x$X1.d2+(1:x$X2.d2)


if(eq==1){ X.int <- x$X1
           X.noi <- x$X2
           ind.int <- ind[[1]]
           ind.noi <- ind[[2]] 
           etap.noi <- x$eta2
                 
}

if(eq==2){ X.int <- x$X2
           X.noi <- x$X1
           ind.int <- ind[[2]]
           ind.noi <- ind[[1]] 
           etap.noi  <- x$eta1
            
}


              
coef.int <- as.numeric(coef(x)[ind.int])
coef.noi <- as.numeric(coef(x)[ind.noi])        

X.int <- X.int[good,]       
X.noi <- X.noi[good,]   
                             
d0 <- d1 <- X.int
d0[,nm.bin] <- 0
d1[,nm.bin] <- 1

eti1 <- d1%*%coef.int 
eti0 <- d0%*%coef.int 
etno <- etap.noi      


if(E==TRUE) ind.excl <- rep(TRUE,x$n) else{ 

      if(treat==TRUE) ind.excl <- as.logical(X.int[,nm.bin]) 
                 else ind.excl <- as.logical(X.int[,nm.bin])!=TRUE
      
    } 



if(x$PL=="P"){
p.int1 <- pmax(pnorm(eti1), 1000*.Machine$double.eps )
p.int0 <- pmax(pnorm(eti0), 1000*.Machine$double.eps )

d.int1 <- dnorm(eti1)
d.int0 <- dnorm(eti0)

p.etn  <- pmax(pnorm(etno), 1000*.Machine$double.eps )
d.etn  <- dnorm(etno)
}




if(x$PL=="PP"){

if(eq==1){

p.int1 <- pmax(pnorm(eti1), 1000*.Machine$double.eps )^x$xi1
p.int0 <- pmax(pnorm(eti0), 1000*.Machine$double.eps )^x$xi1

d.int1 <- pmax(pnorm(eti1), 1000*.Machine$double.eps )^(x$xi1 - 1) * (x$xi1 * dnorm(eti1))
d.int0 <- pmax(pnorm(eti0), 1000*.Machine$double.eps )^(x$xi1 - 1) * (x$xi1 * dnorm(eti0))

p.etn  <- pmax(pnorm(etno), 1000*.Machine$double.eps )^x$xi2   
d.etn  <- pmax(pnorm(etno), 1000*.Machine$double.eps )^(x$xi2 - 1) * (x$xi2 * dnorm(etno))

         }   

if(eq==2){

p.int1 <- pmax(pnorm(eti1), 1000*.Machine$double.eps )^x$xi2
p.int0 <- pmax(pnorm(eti0), 1000*.Machine$double.eps )^x$xi2

d.int1 <- pmax(pnorm(eti1), 1000*.Machine$double.eps )^(x$xi2 - 1) * (x$xi2 * dnorm(eti1))
d.int0 <- pmax(pnorm(eti0), 1000*.Machine$double.eps )^(x$xi2 - 1) * (x$xi2 * dnorm(eti0))

p.etn  <- pmax(pnorm(etno), 1000*.Machine$double.eps )^x$xi1   
d.etn  <- pmax(pnorm(etno), 1000*.Machine$double.eps )^(x$xi1 - 1) * (x$xi1 * dnorm(etno))
        }

}



if(x$PL=="RPP"){

if(eq==1){

p.int1 <- 1-pmax(pnorm(-eti1), 1000*.Machine$double.eps )^x$xi1
p.int0 <- 1-pmax(pnorm(-eti0), 1000*.Machine$double.eps )^x$xi1

d.int1 <- pmax(pnorm(-eti1), 1000*.Machine$double.eps )^(x$xi1 - 1) * (x$xi1 * dnorm(-eti1))
d.int0 <- pmax(pnorm(-eti0), 1000*.Machine$double.eps )^(x$xi1 - 1) * (x$xi1 * dnorm(-eti0))

p.etn  <- 1-pmax(pnorm(-etno), 1000*.Machine$double.eps )^x$xi2   
d.etn  <- pmax(pnorm(-etno), 1000*.Machine$double.eps )^(x$xi2 - 1) * (x$xi2 * dnorm(-etno))

         }   

if(eq==2){

p.int1 <- 1-pmax(pnorm(-eti1), 1000*.Machine$double.eps )^x$xi2
p.int0 <- 1-pmax(pnorm(-eti0), 1000*.Machine$double.eps )^x$xi2

d.int1 <- pmax(pnorm(-eti1), 1000*.Machine$double.eps )^(x$xi2 - 1) * (x$xi2 * dnorm(-eti1))
d.int0 <- pmax(pnorm(-eti0), 1000*.Machine$double.eps )^(x$xi2 - 1) * (x$xi2 * dnorm(-eti0))

p.etn  <- 1-pmax(pnorm(-etno), 1000*.Machine$double.eps )^x$xi1   
d.etn  <- pmax(pnorm(-etno), 1000*.Machine$double.eps )^(x$xi1 - 1) * (x$xi1 * dnorm(-etno))
        }

}


  if(x$PL=="SN"){  
  
  del1 <- -x$xi1/sqrt(1+x$xi1^2)
  del2 <- -x$xi2/sqrt(1+x$xi2^2)


  if(eq==1){

p.int1 <- 2*pmax(pbinorm( eti1, 0, cov12=del1), 1000*.Machine$double.eps )
p.int0 <- 2*pmax(pbinorm( eti0, 0, cov12=del1), 1000*.Machine$double.eps )

d.int1 <- 2*dnorm(eti1)*pmax(pnorm(x$xi1*eti1), 1000*.Machine$double.eps )
d.int0 <- 2*dnorm(eti0)*pmax(pnorm(x$xi1*eti0), 1000*.Machine$double.eps )

p.etn  <- 2*pmax(pbinorm( etno, 0, cov12=del2), 1000*.Machine$double.eps )
d.etn  <- 2*dnorm(etno)*pmax(pnorm(x$xi2*etno), 1000*.Machine$double.eps )

         }   

if(eq==2){

p.int1 <- 2*pmax(pbinorm( eti1, 0, cov12=del2), 1000*.Machine$double.eps )
p.int0 <- 2*pmax(pbinorm( eti0, 0, cov12=del2), 1000*.Machine$double.eps )

d.int1 <- 2*dnorm(eti1)*pmax(pnorm(x$xi2*eti1), 1000*.Machine$double.eps )
d.int0 <- 2*dnorm(eti0)*pmax(pnorm(x$xi2*eti0), 1000*.Machine$double.eps )

p.etn  <- 2*pmax(pbinorm( etno, 0, cov12=del1), 1000*.Machine$double.eps )
d.etn  <- 2*dnorm(etno)*pmax(pnorm(x$xi1*etno), 1000*.Machine$double.eps )
        }    
      
  
   
  }
  

if(x$BivD %in% c("N", "T") ) {ass.p <- x$rho; ass.pst <- coef(x)["athrho"] } else{ ass.p <- x$theta; ass.pst <- coef(x)["theta.star"]} 
 
   if(x$BivD=="N") { 
                     pn.int1 <- qnorm(p.int1)
                     pn.int0 <- qnorm(p.int0)
                     pn.etn  <- qnorm(p.etn) 
                     C.11  <- pmax(pbinorm(pn.int1,pn.etn,cov12=ass.p), 1000*.Machine$double.eps )
                     C.10  <- p.int0 - pmax(pbinorm(pn.int0,pn.etn,cov12=ass.p), 1000*.Machine$double.eps )

                    }else{ 
                    
                     if(x$BivD=="T") pa2 <- x$nu else pa2 <- 0
                     C.11  <- pmax( BiCopCDF(p.int1,p.etn, x$nC, par=ass.p,par2=pa2) , 1000*.Machine$double.eps )
                     C.10  <- p.int0 - pmax( BiCopCDF(p.int0,p.etn, x$nC, par=ass.p, par2=pa2) , 1000*.Machine$double.eps )
                         }           

   est.AT <- mean(   (C.11/p.etn - C.10/(1-p.etn))[ind.excl],na.rm = TRUE   )



            
if(delta==TRUE){

    if(x$BivD %in% c("N","T")      ) add.b <- 1/cosh(coef(x)["athrho"])^2
    if(x$BivD=="F")                  add.b <- 1
    if(x$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) add.b <-  exp(coef(x)["theta.star"])
    if(x$BivD %in% c("C90","C270","J90","J270","G90","G270") ) add.b <- -exp(coef(x)["theta.star"])
   
   dC1 <- copgHs(p1=p.etn,p2=p.int1,teta=ass.p,teta.st=ass.pst,BivD=x$BivD,nC=x$nC,nu=x$nu,xi1=NULL,xi2=NULL,eta1=etno,eta2=eti1,PL=x$PL,eqPL=x$eqPL)
   dC0 <- copgHs(p1=p.etn,p2=p.int0,teta=ass.p,teta.st=ass.pst,BivD=x$BivD,nC=x$nC,nu=x$nu,xi1=NULL,xi2=NULL,eta1=etno,eta2=eti0,PL=x$PL,eqPL=x$eqPL)

   dATT.noint <- ( (dC1$c.copula.be1*p.etn-C.11)/p.etn^2 + (dC0$c.copula.be1*(1-p.etn)-C.10)/(1-p.etn)^2)*d.etn 
   dATT.noint <- colMeans( (c(dATT.noint)*X.noi)[ind.excl,] ) 

   dATT.int   <- colMeans( (c( (dC1$c.copula.be2*d.int1)/p.etn )*d1)[ind.excl,] ) - colMeans( (c( (1-dC0$c.copula.be2)*d.int0/(1-p.etn) )*d0)[ind.excl,])

   dATT.tet   <- mean(((dC1$c.copula.theta/add.b)/p.etn + (dC0$c.copula.theta/add.b)/(1-p.etn))[ind.excl])  


   if(eq==2) dATT <- c(dATT.noint,dATT.int,dATT.tet) else dATT <- c(dATT.int,dATT.noint,dATT.tet) 


   var <- bprobgHs(params=coef(x), PL=x$PL, eqPL=x$eqPL, 
                         respvec=x$respvec, VC=x$VC, 
                         sp=x$sp, qu.mag=x$qu.mag, AT=TRUE)$hessian


   var.eig <- eigen(var, symmetric=TRUE)
    
   if(min(var.eig$values) < .Machine$double.eps){ var <- as.matrix( nearPD( var, ensureSymmetry = FALSE )$mat )
                                                  var.eig <- eigen(var, symmetric=TRUE) }
    
   var <- var.eig$vec%*%tcrossprod(diag(1/var.eig$val),var.eig$vec)  

   delta.AT <- sqrt( t(dATT)%*%var%*%dATT )
                       
}

if(delta==FALSE){


bs <- rmvnorm(n.sim, mean = coef(x), sigma=x$Vb, method=s.meth)


 eti1s <- d1%*%t(bs[,ind.int])    
 eti0s <- d0%*%t(bs[,ind.int])   
 etnos <- X.noi%*%t(bs[,ind.noi]) 
                                  
if(x$PL=="P"){
	p.int1s <- pmax(pnorm(eti1s), 1000*.Machine$double.eps )
	p.int0s <- pmax(pnorm(eti0s), 1000*.Machine$double.eps )
	p.etns  <- pmax(pnorm(etnos), 1000*.Machine$double.eps )
             }


if(x$PL=="PP"){                    
 if(eq==1){
  p.int1s <- pmax(pnorm(eti1s), 1000*.Machine$double.eps )^x$xi1  # this does not account for xi variability, provided it makes sense ... 
  p.int0s <- pmax(pnorm(eti0s), 1000*.Machine$double.eps )^x$xi1
  p.etns  <- pmax(pnorm(etnos), 1000*.Machine$double.eps )^x$xi2 
           }   
 if(eq==2){
  p.int1s <- pmax(pnorm(eti1s), 1000*.Machine$double.eps )^x$xi2
  p.int0s <- pmax(pnorm(eti0s), 1000*.Machine$double.eps )^x$xi2
  p.etns  <- pmax(pnorm(etnos), 1000*.Machine$double.eps )^x$xi1 
           }
}


if(x$PL=="RPP"){
 if(eq==1){
  p.int1s <- 1-pmax(pnorm(-eti1s), 1000*.Machine$double.eps )^x$xi1
  p.int0s <- 1-pmax(pnorm(-eti0s), 1000*.Machine$double.eps )^x$xi1
  p.etns  <- 1-pmax(pnorm(-etnos), 1000*.Machine$double.eps )^x$xi2  
          }   
 if(eq==2){
  p.int1s <- 1-pmax(pnorm(-eti1s), 1000*.Machine$double.eps )^x$xi2
  p.int0s <- 1-pmax(pnorm(-eti0s), 1000*.Machine$double.eps )^x$xi2
  p.etns  <- 1-pmax(pnorm(-etnos), 1000*.Machine$double.eps )^x$xi1  
          }
}


  
if(x$PL=="SN"){  

dim1 <- dim(eti1s)[1]
dim2 <- dim(eti1s)[2]
      
  if(eq==1){

p.int1s <- matrix(2*pmax(pbinorm( as.vector(eti1s), 0, cov12=del1), 1000*.Machine$double.eps ),dim1,dim2)
p.int0s <- matrix(2*pmax(pbinorm( as.vector(eti0s), 0, cov12=del1), 1000*.Machine$double.eps ),dim1,dim2)
p.etns  <- matrix(2*pmax(pbinorm( as.vector(etnos), 0, cov12=del2), 1000*.Machine$double.eps ),dim1,dim2)

         }   

if(eq==2){

p.int1s <- matrix(2*pmax(pbinorm( as.vector(eti1s), 0, cov12=del2), 1000*.Machine$double.eps ),dim1,dim2)
p.int0s <- matrix(2*pmax(pbinorm( as.vector(eti0s), 0, cov12=del2), 1000*.Machine$double.eps ),dim1,dim2)
p.etns  <- matrix(2*pmax(pbinorm( as.vector(etnos), 0, cov12=del1), 1000*.Machine$double.eps ),dim1,dim2)

        }    
   
  }  
  
  


if(x$BivD!="N"){

   if(x$BivD=="F")                   ass.ps <- bs[,p.rho] + epsilon
   if(x$BivD=="T")                  {ass.ps <- tanh(bs[,p.rho]); ass.ps <- ifelse(ass.ps %in% c(-1,1), sign(ass.ps)*0.9999999, ass.ps) }
  
   if(x$BivD %in% c("C0", "C180") )  ass.ps <-   exp(bs[,p.rho]) + epsilon
   if(x$BivD %in% c("C90","C270") )  ass.ps <- -(exp(bs[,p.rho]) + epsilon)

   if(x$BivD %in% c("J0", "J180") )  ass.ps <-    1 + exp(bs[,p.rho]) + epsilon
   if(x$BivD %in% c("J90","J270") )  ass.ps <- -( 1 + exp(bs[,p.rho]) + epsilon)
 
   if(x$BivD %in% c("G0", "G180") )  ass.ps <-    1 + exp(bs[,p.rho])
   if(x$BivD %in% c("G90","G270") )  ass.ps <- -( 1 + exp(bs[,p.rho]) )
   
   ass.ps <- ifelse(ass.ps=="Inf" ,  2.688117e+43, ass.ps )
   ass.ps <- ifelse(ass.ps=="-Inf", -2.688117e+43, ass.ps )

}

if(x$BivD=="T") asp2 <- rep(x$nu,n.sim) else asp2 <- rep(0,n.sim) 





for(i in 1:n.sim){ 

 if(x$BivD=="N"){    pn.int1s <- qnorm(p.int1s[,i])
                     pn.int0s <- qnorm(p.int0s[,i])
                     pn.etns  <- qnorm(p.etns[,i]) 

est.ATb[i] <- mean(   (pmax(pbinorm(pn.int1s,pn.etns,cov12=tanh(bs[i,p.rho])), 1000*.Machine$double.eps)/p.etns[,i] - (p.int0s[,i] - pmax(pbinorm(pn.int0s,pn.etns,cov12=tanh(bs[i,p.rho])), 1000*.Machine$double.eps))/(1-p.etns[,i]))[ind.excl],na.rm = TRUE   )

                }else{

 C.11 <- pmax( BiCopCDF(p.int1s[,i],p.etns[,i], x$nC, par=ass.ps[i],par2=asp2[i]) , 1000*.Machine$double.eps )
 C.10 <- p.int0s[,i] - pmax( BiCopCDF(p.int0s[,i],p.etns[,i], x$nC, par=ass.ps[i],par2=asp2[i]) , 1000*.Machine$double.eps )

 est.ATb[i] <- mean(   (C.11/p.etns[,i] - C.10/(1-p.etns[,i]))[ind.excl],na.rm = TRUE   )
                     }


                 }



}  # end big loop   



            

if(delta==TRUE){esd.AT <- delta.AT*qnorm(prob.lev/2,lower.tail = FALSE) 
                   CIs <- c(est.AT - esd.AT, est.AT + esd.AT)
               }else CIs <- as.numeric(quantile(est.ATb,c(prob.lev/2,1-prob.lev/2),na.rm = TRUE))

res <- c(CIs[1],est.AT,CIs[2])
out <- list(res=res, prob.lev=prob.lev, est.ATb=est.ATb)
 
class(out) <- "AT"

out

}




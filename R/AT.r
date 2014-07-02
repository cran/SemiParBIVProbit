AT <- function(x,eq,nm.bin="",E=TRUE,treat=TRUE,delta=FALSE, prob.lev=0.05,s.meth="svd",n.sim=1000){


tpc <- c("BB1.0","BB1.180","BB1.90","BB1.270",
         "BB6.0","BB6.180","BB6.90","BB6.270",
         "BB7.0","BB7.180","BB7.90","BB7.270",
         "BB8.0","BB8.180","BB8.90","BB8.270")


if(is.character(nm.bin)==FALSE) stop("nm.bin is not a character!")
if(x$sel==TRUE) stop("Calculation of average treatment effects for sample selection models not implemented.")
if(x$PL!="P" || x$BivD %in% tpc ) delta <- FALSE

est.ATb <- NA
ass.ps2 <- ass.p2 <- NULL
ind <- list()

p.rho <- length(coef(x))
if(x$BivD %in% tpc ) p.rho <- p.rho - 1

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
  

if(x$BivD %in% c("N", "T") ) {ass.p <- x$rho; ass.pst <- coef(x)["athrho"] } else{ ass.p <- x$theta; ass.pst <- coef(x)["theta.star"]; 
                                                                                   ass.p2 <- x$delta} 
 
   if(x$BivD=="N") { 
                     pn.int1 <- qnorm(p.int1)
                     pn.int0 <- qnorm(p.int0)
                     pn.etn  <- qnorm(p.etn) 
                     C.11  <- abs(pbinorm(pn.int1,pn.etn,cov12=ass.p))
                     C.10  <- p.int0 - abs(pbinorm(pn.int0,pn.etn,cov12=ass.p))

                    }else{ 
                    
                     if(x$BivD=="T") pa2 <- x$nu else {if(!is.null(ass.p2)) pa2 <- ass.p2 else pa2 <- 0} 
                     C.11  <- pmax( BiCopCDF(p.int1,p.etn, x$nC, par=ass.p,par2=pa2) , 1000*.Machine$double.eps )
                     C.10  <- p.int0 - pmax( BiCopCDF(p.int0,p.etn, x$nC, par=ass.p, par2=pa2) , 1000*.Machine$double.eps )
                         }           

   est.AT <- mean(   (C.11/p.etn - C.10/(1-p.etn))[ind.excl],na.rm = TRUE   )



            
if(delta==TRUE){

    if(x$BivD %in% c("N","T")      ) add.b <- 1/cosh(coef(x)["athrho"])^2
    if(x$BivD=="F")                  add.b <- 1
    if(x$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) add.b <-  exp(coef(x)["theta.star"])
    if(x$BivD %in% c("C90","C270","J90","J270","G90","G270") ) add.b <- -exp(coef(x)["theta.star"])


   dC1 <- copgHs(p1=p.etn,p2=p.int1,teta=ass.p,teta.st=ass.pst,BivD=x$BivD,nC=x$nC,nu=x$nu,xi1=x$xi1,xi1.st=log(x$xi1),xi2=x$xi2,xi2.st=log(x$xi2),PL=x$PL,eta1=etno,eta2=eti1,eqPL=x$eqPL)
   dC0 <- copgHs(p1=p.etn,p2=p.int0,teta=ass.p,teta.st=ass.pst,BivD=x$BivD,nC=x$nC,nu=x$nu,xi1=x$xi1,xi1.st=log(x$xi1),xi2=x$xi2,xi2.st=log(x$xi2),PL=x$PL,eta1=etno,eta2=eti0,eqPL=x$eqPL)

   dATT.noint <- ( (dC1$c.copula.be1*p.etn-C.11)/p.etn^2 + (dC0$c.copula.be1*(1-p.etn)-C.10)/(1-p.etn)^2)*d.etn 
   dATT.noint <- colMeans( (c(dATT.noint)*X.noi)[ind.excl,] ) 

   dATT.int   <- colMeans( (c( (dC1$c.copula.be2*d.int1)/p.etn )*d1)[ind.excl,] ) - colMeans( (c( (1-dC0$c.copula.be2)*d.int0/(1-p.etn) )*d0)[ind.excl,])

   dATT.tet   <- mean(((dC1$c.copula.theta/add.b)/p.etn + (dC0$c.copula.theta/add.b)/(1-p.etn))[ind.excl])  


   if(eq==2) dATT <- c(dATT.noint,dATT.int,dATT.tet) else dATT <- c(dATT.int,dATT.noint,dATT.tet) 

#if(x$PL!="P"){
#
#  add.bl1 <- x$lambda1 
#   add.bl2 <- x$lambda2
#
#   dATT.l1    <- mean((  ((dC1$c.copula.lambda1/add.bl1)*p.etn - C.11*(dC1$der.p1.lambda1/add.bl1))/p.etn^2 - 
#                         ((-dC0$c.copula.lambda1/add.bl1)*(1-p.etn) + C.10*(dC0$der.p1.lambda1/add.bl1))/(1-p.etn)^2)[ind.excl]) 
#                  
#   dATT.l2    <- mean((  (dC1$c.copula.lambda2/add.bl2)/p.etn - ((dC0$der.p2.lambda2/add.bl2)-(dC0$c.copula.lambda2/add.bl2))/(1-p.etn) )[ind.excl]) 
#
#}
#  if(x$PL=="P"){  
#               }
#
#  if(x$PL!="P"){  
#
#       if(eq==2) dATT <- c(dATT.noint,dATT.int,dATT.tet,dATT.l1,dATT.l2) else dATT <- c(dATT.int,dATT.noint,dATT.tet,dATT.l2,dATT.l1)
#
#       l.da <- length(dATT); if(x$eqPL=="first") l.da <- l.da - 1 else l.da <- l.da  
#       if(x$eqPL!="both") dATT <- dATT[-c(l.da)] 
#
#               }
#var <- x$Vb   
#if(x$PL=="P") fun.AT <- bprobgHs else fun.AT <- bprobgHsPL   


   var <- solve(bprobgHs(params=coef(x), BivD=x$BivD, nC=x$nC, nu=x$nu, PL=x$PL, eqPL=x$eqPL, H.n=TRUE, 
            y1.y2=x$y1.y2, y1.cy2=x$y1.cy2, cy1.y2=x$cy1.y2, cy1.cy2=x$cy1.cy2, cy1=x$cy1, 
            X1=x$X1, X2=x$X2, weights=x$weights, X1.d2=x$X1.d2, X2.d2=x$X2.d2, pPen1=x$pPen1, pPen2=x$pPen2, 
            sp=x$sp, qu.mag=x$qu.mag, gp1=x$gp1, gp2=x$gp2, fp=x$fp, l.sp1=x$l.sp1, l.sp2=x$l.sp2, AT=TRUE)$hessian)



   #if(x$PL!="P" && x$eqPL!="both") var <- adiag(var,0)   
   #   ; if(x$PL!="P") var <- var[1:p.rho,1:p.rho]    # x$Vb%*%x$HeSh%*%x$Vb


   delta.AT <- sqrt( t(dATT)%*%var%*%dATT )
                       
}

if(delta==FALSE){



bs <- rmvnorm(n.sim, mean = coef(x), sigma=x$Vb, method=s.meth)


 eti1s <- d1%*%t(bs[,ind.int])    
 eti0s <- d0%*%t(bs[,ind.int])   
 etnos <- X.noi%*%t(bs[,ind.noi]) 
                                  
if(x$PL=="P"){
	p.int1s <- pnorm(eti1s)
	p.int0s <- pnorm(eti0s)
	p.etns  <- pnorm(etnos)
             }


#if(x$PL!="P"){
#  if(x$eqPL=="both"){
#                     lambda1s <- exp(bs[i,(p.rho+1)]) + epsilon
#                     lambda2s <- exp(bs[i,(p.rho+2)]) + epsilon
#                     }
#  if(x$eqPL=="first"){
#                     lambda1s <- exp(bs[i,(p.rho+1)]) + epsilon
#                     lambda2s <- 1
#                     }
#  if(x$eqPL=="second"){
#                     lambda1s <- 1
#                     lambda2s <- exp(bs[i,(p.rho+1)]) + epsilon
#                     } 
#}


if(x$PL=="PP"){                    
 if(eq==1){
  p.int1s <- pnorm(eti1s)^x$xi1  # here I need to replace with bit above to account for xi variability (provided it makes sense ... )
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



#if(x$PL=="SN"){  
#      
#dim1 <- dim(eti1s)[1]
#dim2 <- dim(eti1s)[2]
#      
#      
#  if(eq==1){
#
#p.int1s <- matrix(psn(eti1s, alpha=x$xi1),dim1,dim2)
#p.int0s <- matrix(psn(eti0s, alpha=x$xi1),dim1,dim2)
#p.etns  <- matrix(psn(etnos, alpha=x$xi2),dim1,dim2)
#
#         }   
#
#if(eq==2){
#
#p.int1s <- matrix(psn(eti1s, alpha=x$xi2),dim1,dim2)
#p.int0s <- matrix(psn(eti0s, alpha=x$xi2),dim1,dim2)
#p.etns  <- matrix(psn(etnos, alpha=x$xi1),dim1,dim2)
#
#        }    
#   
#  }
  
  
  
  
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
  
  
  
  
  
  
  
  
  


  #p.int1s <- ifelse(p.int1s==0,0.000000001,p.int1s)
  #p.int1s <- ifelse(p.int1s==1,0.999999999,p.int1s)
  #p.int0s <- ifelse(p.int0s==0,0.000000001,p.int0s)
  #p.int0s <- ifelse(p.int0s==1,0.999999999,p.int0s)
  #p.etns  <- ifelse(p.etns==1,0.999999999,p.etns)
  #p.etns  <- ifelse(p.etns==1,0.999999999,p.etns)


if(x$BivD!="N"){

   if(x$BivD=="F")                   ass.ps <- bs[,p.rho] + epsilon
   if(x$BivD=="T")                  {ass.ps <- tanh(bs[,p.rho]); ass.ps <- ifelse(ass.ps %in% c(-1,1), sign(ass.ps)*0.9999999, ass.ps) }
  
   if(x$BivD %in% c("C0", "C180") )  ass.ps <-   exp(bs[,p.rho]) + epsilon
   if(x$BivD %in% c("C90","C270") )  ass.ps <- -(exp(bs[,p.rho]) + epsilon)

   if(x$BivD %in% c("J0", "J180") )  ass.ps <-    1 + exp(bs[,p.rho]) + epsilon
   if(x$BivD %in% c("J90","J270") )  ass.ps <- -( 1 + exp(bs[,p.rho]) + epsilon)
 
   if(x$BivD %in% c("G0", "G180") )  ass.ps <-    1 + exp(bs[,p.rho])
   if(x$BivD %in% c("G90","G270") )  ass.ps <- -( 1 + exp(bs[,p.rho]) )
   
   if(x$BivD %in% c("BB1.0", "BB1.180")){ass.ps <-   exp(bs[,p.rho]) + epsilon;  ass.ps2 <-   exp(bs[,p.rho+1]) + 1}
   if(x$BivD %in% c("BB1.90","BB1.270")){ass.ps <- -(exp(bs[,p.rho]) + epsilon); ass.ps2 <- -(exp(bs[,p.rho+1]) + 1)}

   if(x$BivD %in% c("BB6.0", "BB6.180")){ass.ps <-   exp(bs[,p.rho]) + 1; ass.ps2 <-   exp(bs[,p.rho+1]) + 1}
   if(x$BivD %in% c("BB6.90","BB6.270")){ass.ps <- -(exp(bs[,p.rho]) + 1);ass.ps2 <- -(exp(bs[,p.rho+1]) + 1)}

   if(x$BivD %in% c("BB7.0", "BB7.180")){ass.ps <-   exp(bs[,p.rho]) + 1; ass.ps2 <-   exp(bs[,p.rho+1]) + epsilon}
   if(x$BivD %in% c("BB7.90","BB7.270")){ass.ps <- -(exp(bs[,p.rho]) + 1);ass.ps2 <- -(exp(bs[,p.rho+1]) + epsilon)}

   if(x$BivD %in% c("BB8.0", "BB8.180")){ass.ps <-   exp(bs[,p.rho]) + 1; ass.ps2 <-  pnorm(bs[,p.rho+1]); ass.ps2 <- ifelse(ass.ps2==0,  epsilon, ass.ps2)}
   if(x$BivD %in% c("BB8.90","BB8.270")){ass.ps <- -(exp(bs[,p.rho]) + 1);ass.ps2 <- -pnorm(bs[,p.rho+1]); ass.ps2 <- ifelse(ass.ps2==0, -epsilon, ass.ps2)}   
   
   ass.ps <- ifelse(ass.ps=="Inf" ,  2.688117e+43, ass.ps )
   ass.ps <- ifelse(ass.ps=="-Inf", -2.688117e+43, ass.ps )


   #if(ass.ps2=="Inf" || ass.ps2=="-Inf" ) ass.ps2 <- sign(ass.ps2)*2.688117e+43

}

if(x$BivD=="T") asp2 <- rep(x$nu,n.sim) else {if(is.null(ass.ps2)) asp2 <- rep(0,n.sim) else asp2 <- ass.ps2} 





for(i in 1:n.sim){ 

 if(x$BivD=="N"){    pn.int1s <- qnorm(p.int1s[,i])
                     pn.int0s <- qnorm(p.int0s[,i])
                     pn.etns  <- qnorm(p.etns[,i]) 

est.ATb[i] <- mean(   (abs(pbinorm(pn.int1s,pn.etns,cov12=tanh(bs[i,p.rho])))/p.etns[,i] - (p.int0s[,i] - abs(pbinorm(pn.int0s,pn.etns,cov12=tanh(bs[i,p.rho]))))/(1-p.etns[,i]))[ind.excl],na.rm = TRUE   )

                }else{

 C.11 <- pmax( BiCopCDF(p.int1s[,i],p.etns[,i], x$nC, par=ass.ps[i],par2=asp2[i]) , 1000*.Machine$double.eps )
 C.10 <- p.int0s[,i] - pmax( BiCopCDF(p.int0s[,i],p.etns[,i], x$nC, par=ass.ps[i],par2=asp2[i]) , 1000*.Machine$double.eps )

 est.ATb[i] <- mean(   (C.11/p.etns[,i] - C.10/(1-p.etns[,i]))[ind.excl],na.rm = TRUE   )
                     }


                 }



}  # end general loop   



            

if(delta==TRUE){esd.AT <- delta.AT*qnorm(prob.lev/2,lower.tail = FALSE) 
                   CIs <- c(est.AT - esd.AT, est.AT + esd.AT)
               }else CIs <- as.numeric(quantile(est.ATb,c(prob.lev/2,1-prob.lev/2),na.rm = TRUE))

res <- c(CIs[1],est.AT,CIs[2])
out <- list(res=res, prob.lev=prob.lev, est.ATb=est.ATb)
 
class(out) <- "AT"

out

}




bprobgHsContUniv3 <- function(params, respvec, VC, ps, AT = FALSE){

epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)

weights <- VC$weights

eta2 <- VC$X1%*%params[1:VC$X1.d2] # this is eta1
eta2 <- eta.tr(eta2, VC$margins[1])
    
if(is.null(VC$X2))  sigma2.st <- params[(VC$X1.d2 + 1)] 
if(!is.null(VC$X2)) sigma2.st <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]

if(is.null(VC$X3))  nu.st <- params[(VC$X1.d2 + 2)] 
if(!is.null(VC$X3)) nu.st <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]

sstr1 <- esp.tr(sigma2.st, VC$margins[1])  
sigma2.st <- sstr1$vrb.st 
sigma2    <- sstr1$vrb 
    
sstr1 <- enu.tr(nu.st, VC$margins[1])  
nu.st <- sstr1$vrb.st 
nu    <- sstr1$vrb           
 
if(VC$margins[1] %in% VC$m3)  dHs <-      distrHs(respvec$y1, eta2, sigma2, sigma2.st, nu, nu.st, margin2=VC$margins[1], naive = TRUE)
if(VC$margins[1] %in% VC$m3d) dHs <- distrHsDiscr(respvec$y1, eta2, sigma2, sigma2.st, nu, nu.st, margin2=VC$margins[2], naive = TRUE, y2m = VC$y1m)


########################################################################################################

pdf2                         <- dHs$pdf2
derpdf2.dereta2              <- dHs$derpdf2.dereta2 
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st  
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2         
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  


der2pdf2.dereta2dernu.st    <- dHs$der2pdf2.dereta2dernu.st   
der2pdf2.dersigma2.stdernu.st <- dHs$der2pdf2.sigma2.st2dernu.st
derpdf2.dernu.st            <- dHs$derpdf2.dernu.st           
der2pdf2.dernu.st2          <- dHs$der2pdf2.dernu.st2         
    
########################################################################################################

if(VC$robust == FALSE){

l.par <- weights*log(pdf2)
res   <- -sum(l.par)
d.psi <- 1

}else{

l.par1    <- log(pdf2)
Robj.lpar <- llpsi(l.par1, VC$rc)
psi       <- Robj.lpar$psi
d.psi     <- Robj.lpar$d.psi
d2.psi    <- Robj.lpar$d2.psi 
l.par     <- weights*( psi )  # weight in or out??
res       <- -sum(l.par)
}

########################################################################################################
 
dl.dbe0       <- derpdf2.dereta2/pdf2 
dl.dsigma.st0 <- derpdf2.dersigma2.st/pdf2 
dl.dnu.st0    <- derpdf2.dernu.st/pdf2   
 
dl.dbe       <- -weights*d.psi*dl.dbe0      
dl.dsigma.st <- -weights*d.psi*dl.dsigma.st0
dl.dnu.st    <- -weights*d.psi*dl.dnu.st0   
  
  
if(VC$robust == FALSE){   
   
d2l.be.be       <- -weights*( (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 )
d2l.sigma.sigma <- -weights*( (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 )
d2l.be.sigma    <- -weights*( (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 )
d2l.be.nu       <- -weights*( (der2pdf2.dereta2dernu.st*pdf2 - derpdf2.dereta2*derpdf2.dernu.st)/pdf2^2 )
d2l.sigma.nu    <- -weights*( (der2pdf2.dersigma2.stdernu.st*pdf2-(derpdf2.dersigma2.st*derpdf2.dernu.st))/pdf2^2 )
d2l.nu.nu       <- -weights*( (der2pdf2.dernu.st2*pdf2-(derpdf2.dernu.st)^2)/pdf2^2 )
  
}  



if(VC$robust == TRUE){   
   
d2l.be.be0       <- (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 
d2l.sigma.sigma0 <- (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 
d2l.be.sigma0    <- (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 
d2l.be.nu0       <- (der2pdf2.dereta2dernu.st*pdf2 - derpdf2.dereta2*derpdf2.dernu.st)/pdf2^2 
d2l.sigma.nu0    <- (der2pdf2.dersigma2.stdernu.st*pdf2-(derpdf2.dersigma2.st*derpdf2.dernu.st))/pdf2^2 
d2l.nu.nu0       <- (der2pdf2.dernu.st2*pdf2-(derpdf2.dernu.st)^2)/pdf2^2 
  
d2l.be.be       <- -weights*(d2.psi*dl.dbe0^2                + d.psi*d2l.be.be0          )
d2l.sigma.sigma <- -weights*(d2.psi*dl.dsigma.st0^2          + d.psi*d2l.sigma.sigma0    )
d2l.be.sigma    <- -weights*(d2.psi*dl.dbe0*dl.dsigma.st0    + d.psi*d2l.be.sigma0       ) 
d2l.be.nu       <- -weights*(d2.psi*dl.dbe0*dl.dnu.st0       + d.psi*d2l.be.nu0          )
d2l.sigma.nu    <- -weights*(d2.psi*dl.dsigma.st0*dl.dnu.st0 + d.psi*d2l.sigma.nu0       )
d2l.nu.nu       <- -weights*(d2.psi*dl.dnu.st0^2             + d.psi*d2l.nu.nu0          )
   
}  



     
     
########################################################################################################
     


if( is.null(VC$X2) ){

  G   <- c( colSums( c(dl.dbe)*VC$X1 ) ,
            sum( dl.dsigma.st ), sum( dl.dnu.st ) )
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be.sigma))))) 
  be.nu    <- t(t(rowSums(t(VC$X1*c(d2l.be.nu))))) 

  H <- rbind( cbind( be.be      , be.sigma,             be.nu             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma), sum(d2l.sigma.nu) ),
              cbind( t(be.nu),    sum(d2l.sigma.nu),    sum(d2l.nu.nu)    )  )     
              
              }


if( !is.null(VC$X2) ){

  G   <- c( colSums( c(dl.dbe)*VC$X1      ),
            colSums( c(dl.dsigma.st)*VC$X2),
            colSums( c(dl.dnu.st)*VC$X3   )  )
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- crossprod(VC$X1*c(d2l.be.sigma),VC$X2)
  be.nu    <- crossprod(VC$X1*c(d2l.be.nu),VC$X3)
  
  si.sigma <- crossprod(VC$X2*c(d2l.sigma.sigma),VC$X2)  
  si.nu    <- crossprod(VC$X3*c(d2l.nu.nu),VC$X3)    
  
  sa.nu    <- crossprod(VC$X2*c(d2l.sigma.nu),VC$X3)    
  
  H <- rbind( cbind( be.be      , be.sigma , be.nu ), 
              cbind( t(be.sigma), si.sigma,  sa.nu ),
              cbind( t(be.nu),    t(sa.nu),  si.nu ) )       
              
              }







if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
  S.h  <- ps$S.h  


  if( length(S.h) != 1){
  
  S.h1 <- 0.5*crossprod(params,S.h)%*%params
  S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h   
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2) 
  

list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.par=l.par, 
     ps = ps, sigma2.st = sigma2.st, nu.st = nu.st, etas1 = sigma2.st, etan1 = nu.st, 
     BivD=VC$BivD, eta1 = eta2, eta2 = eta2, sigma2 = sigma2, nu = nu)      

}

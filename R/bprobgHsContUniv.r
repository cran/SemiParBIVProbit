bprobgHsContUniv <- function(params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){

X2 <- X3 <- 1

eta2 <- VC$X2%*%params[1:VC$X2.d2]

epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)

if(is.null(VC$X3))  sigma2.st <- params[(VC$X2.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- VC$X3%*%params[(VC$X2.d2+1):(VC$X2.d2+VC$X3.d2)]

sigma2 <- exp(sigma2.st) + epsilon

y2 <- respvec$y2 
weights <- VC$weights
  
dHs  <- distrHs(y2, eta2, sigma2, sigma2.st, margin2=VC$margins[2], naive = TRUE)
  
pdf2                         <- dHs$pdf2
derpdf2.dereta2              <- dHs$derpdf2.dereta2 
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st  
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2         
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  
  
    
########################################################################################################

l.par <- weights*log(pdf2)
res   <- -sum(l.par)
  
########################################################################################################
 
dl.dbe       <- -weights*( derpdf2.dereta2/pdf2 )
dl.dsigma.st <- -weights*( derpdf2.dersigma2.st/pdf2 )
                     
d2l.be.be        <- -weights*( (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 )
d2l.sigma.sigma  <- -weights*( (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 )
d2l.be.sigma     <- -weights*( (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 )
     
########################################################################################################
     

if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2 ) ,
            sum( dl.dsigma.st )
          )
                

  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be.sigma))))) 

  H <- rbind( cbind( be.be      , be.sigma             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma) )  
            ) 
  
}




if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2        ) ,
            colSums( c(dl.dsigma.st)*VC$X3  )
          )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- crossprod(VC$X2*c(d2l.be.sigma),VC$X3)
  si.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)  

  H <- rbind( cbind( be.be      , be.sigma  ), 
              cbind( t(be.sigma), si.sigma  )  
            )  
   
    
}



if( (VC$l.sp2==0 && VC$l.sp3==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, eq1 = "no")


if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2) 
  
rm(X2, X3)  
  

         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps, sigma2.st = sigma2.st,
              eta2=eta2, 
dl.dbe          = dl.dbe,            
dl.dsigma.st    = dl.dsigma.st,
d2l.be.be       = d2l.be.be,
d2l.sigma.sigma = d2l.sigma.sigma,
d2l.be.sigma    = d2l.be.sigma,
              BivD=VC$BivD)      

}




     
























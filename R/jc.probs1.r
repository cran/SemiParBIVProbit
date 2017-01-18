jc.probs1 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link){

######################################################################################################

epsilon <- 0.0000001 
nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- dof <- p12s <- dofs <- C1s <- C2s <- C11s <- C01s <- C10s <- C00s <- NULL

######################################################################################################


if(type == "bivariate"){ 



if(!missing(newdata)){ #################

nu1 <- nu2 <- dof <- sigma21 <- sigma22 <- NA

eta1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata)
eta2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata)

if( !is.null(x$X3) ){ # X3


if(x$margins[1] %in% c(cont2par, cont3par))                                           sigma21 <- esp.tr(predict.SemiParBIVProbit(x, eq = 3, newdata = newdata), x$margins[1])$vrb
if(x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par, cont3par)) sigma22 <- esp.tr(predict.SemiParBIVProbit(x, eq = 4, newdata = newdata), x$margins[2])$vrb
if(x$margins[1] %in% cont1par && x$margins[2] %in% c(cont2par,cont3par))              sigma22 <- esp.tr(predict.SemiParBIVProbit(x, eq = 3, newdata = newdata), x$margins[2])$vrb


if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3) ){

if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){ eq.nu1 <- 5;  eq.nu2 <- 6;  eq.dof <- 7; eq.th <- 8}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.dof <- 5; eq.th <- 6}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 5;  eq.dof <- 6; eq.th <- 7}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){ eq.nu1 <- 5;  eq.nu2 <- NA; eq.dof <- 6; eq.th <- 7}

dof <- dof.tr(predict.SemiParBIVProbit(x, eq = eq.dof, newdata = newdata))$vao 


                   }else{
                   
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){ eq.nu1 <- 5;  eq.nu2 <- 6;  eq.th <- 7}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.th <- 5}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 5;  eq.th <- 6}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){ eq.nu1 <- 5;  eq.nu2 <- NA; eq.th <- 6} 

if(x$margins[1] %in% cont1par && x$margins[2] %in% cont1par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.th <- 3}
if(x$margins[1] %in% cont1par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.th <- 4}
if(x$margins[1] %in% cont1par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 4 ; eq.th <- 5}
                   
}                   
                   


if(x$margins[1] %in% cont3par) nu1 <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu1, newdata = newdata), x$margins[1])$vrb
if(x$margins[2] %in% cont3par) nu2 <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu2, newdata = newdata), x$margins[2])$vrb

theta <- teta.tr(x$VC, predict.SemiParBIVProbit(x, eq = eq.th, newdata = newdata))$teta

} # X3


if( is.null(x$X3) ){

sigma21 <- x$sigma21
sigma22 <- x$sigma22

nu1 <- x$nu1 
nu2 <- x$nu2

theta <- x$theta 
dof <- x$dof

                   }


} ############



if(missing(newdata)){

eta1 <- x$eta1
eta2 <- x$eta2

sigma21 <- x$sigma21
sigma22 <- x$sigma22

nu1 <- x$nu1 
nu2 <- x$nu2

theta <- x$theta 
dof   <- x$dof

}


if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p1 <- distrHsAT(y1, eta1, sigma21, nu1, x$margins[1])$p2 
if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p2 <- distrHsAT(y2, eta2, sigma22, nu2, x$margins[2])$p2 

if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d)) p1pdf1 <- distrHsATDiscr(y1, eta1, sigma21, nu = 1, x$margins[1], x$VC$y1m)
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)) p2pdf2 <- distrHsATDiscr(y2, eta2, sigma22, nu = 1, x$margins[2], x$VC$y2m)








if(x$BivD %in% x$BivD2){

nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]

p12 <- C11 <- C01 <- C10 <- C00 <- NA
 
 
if(x$margins[1] %in% c(x$VC$m2, x$VC$m3) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){

if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- BiCDF(p1[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof)
if(length(theta) == 1) p12[x$teta.ind1] <- BiCDF(p1[x$teta.ind1], p2[x$teta.ind1], nC1, theta, dof)
                          }                       
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- BiCDF(p1[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof)
if(length(theta) == 1) p12[x$teta.ind2] <- BiCDF(p1[x$teta.ind2], p2[x$teta.ind2], nC2, theta, dof)
                          }                            
                                                                                   }


if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){

if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- BiCDF(p1pdf1$p2[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof) - BiCDF(mm(p1pdf1$p2[x$teta.ind1] - p1pdf1$pdf2[x$teta.ind1]), p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof); p12[x$teta.ind1] <- ifelse(p12[x$teta.ind1] < epsilon, epsilon, p12[x$teta.ind1]) 
if(length(theta) == 1) p12[x$teta.ind1] <- BiCDF(p1pdf1$p2[x$teta.ind1], p2[x$teta.ind1], nC1, theta, dof)              - BiCDF(mm(p1pdf1$p2[x$teta.ind1] - p1pdf1$pdf2[x$teta.ind1]), p2[x$teta.ind1], nC1, theta,              dof); p12[x$teta.ind1] <- ifelse(p12[x$teta.ind1] < epsilon, epsilon, p12[x$teta.ind1]) 
                          }                      
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- BiCDF(p1pdf1$p2[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof) - BiCDF(mm(p1pdf1$p2[x$teta.ind2] - p1pdf1$pdf2[x$teta.ind2]), p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof); p12[x$teta.ind2] <- ifelse(p12[x$teta.ind2] < epsilon, epsilon, p12[x$teta.ind2]) 
if(length(theta) == 1) p12[x$teta.ind2] <- BiCDF(p1pdf1$p2[x$teta.ind2], p2[x$teta.ind2], nC2, theta, dof)              - BiCDF(mm(p1pdf1$p2[x$teta.ind2] - p1pdf1$pdf2[x$teta.ind2]), p2[x$teta.ind2], nC2, theta,              dof); p12[x$teta.ind2] <- ifelse(p12[x$teta.ind2] < epsilon, epsilon, p12[x$teta.ind2]) 
                          } 
                                                                                    }




if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){


if( length(x$teta1) != 0){

if(length(theta) > 1){  

  C11[x$teta.ind1] <- BiCDF(p1pdf1$p2[x$teta.ind1],                              p2pdf2$p2[x$teta.ind1],                              nC1, theta[x$teta.ind1], dof)
  C01[x$teta.ind1] <- BiCDF(mm(p1pdf1$p2[x$teta.ind1]-p1pdf1$pdf2[x$teta.ind1]), p2pdf2$p2[x$teta.ind1],                              nC1, theta[x$teta.ind1], dof)
  C10[x$teta.ind1] <- BiCDF(p1pdf1$p2[x$teta.ind1],                              mm(p2pdf2$p2[x$teta.ind1]-p2pdf2$pdf2[x$teta.ind1]), nC1, theta[x$teta.ind1], dof)
  C00[x$teta.ind1] <- BiCDF(mm(p1pdf1$p2[x$teta.ind1]-p1pdf1$pdf2[x$teta.ind1]), mm(p2pdf2$p2[x$teta.ind1]-p2pdf2$pdf2[x$teta.ind1]), nC1, theta[x$teta.ind1], dof)

  p12[x$teta.ind1] <- C11[x$teta.ind1] - C01[x$teta.ind1] - C10[x$teta.ind1] + C00[x$teta.ind1]
  p12[x$teta.ind1] <- ifelse(p12[x$teta.ind1] < epsilon, epsilon, p12[x$teta.ind1])

                     }


if(length(theta) == 1){  

  C11[x$teta.ind1] <- BiCDF(p1pdf1$p2[x$teta.ind1],                              p2pdf2$p2[x$teta.ind1],                              nC1, theta, dof)
  C01[x$teta.ind1] <- BiCDF(mm(p1pdf1$p2[x$teta.ind1]-p1pdf1$pdf2[x$teta.ind1]), p2pdf2$p2[x$teta.ind1],                              nC1, theta, dof)
  C10[x$teta.ind1] <- BiCDF(p1pdf1$p2[x$teta.ind1],                              mm(p2pdf2$p2[x$teta.ind1]-p2pdf2$pdf2[x$teta.ind1]), nC1, theta, dof)
  C00[x$teta.ind1] <- BiCDF(mm(p1pdf1$p2[x$teta.ind1]-p1pdf1$pdf2[x$teta.ind1]), mm(p2pdf2$p2[x$teta.ind1]-p2pdf2$pdf2[x$teta.ind1]), nC1, theta, dof)

  p12[x$teta.ind1] <- C11[x$teta.ind1] - C01[x$teta.ind1] - C10[x$teta.ind1] + C00[x$teta.ind1]
  p12[x$teta.ind1] <- ifelse(p12[x$teta.ind1] < epsilon, epsilon, p12[x$teta.ind1])

                     } 

                          }
                          
                                                             
if( length(x$teta2) != 0){

if(length(theta) > 1){  

  C11[x$teta.ind2] <- BiCDF(p1pdf1$p2[x$teta.ind2],                              p2pdf2$p2[x$teta.ind2],                              nC2, theta[x$teta.ind2], dof)
  C01[x$teta.ind2] <- BiCDF(mm(p1pdf1$p2[x$teta.ind2]-p1pdf1$pdf2[x$teta.ind2]), p2pdf2$p2[x$teta.ind2],                              nC2, theta[x$teta.ind2], dof)
  C10[x$teta.ind2] <- BiCDF(p1pdf1$p2[x$teta.ind2],                              mm(p2pdf2$p2[x$teta.ind2]-p2pdf2$pdf2[x$teta.ind2]), nC2, theta[x$teta.ind2], dof)
  C00[x$teta.ind2] <- BiCDF(mm(p1pdf1$p2[x$teta.ind2]-p1pdf1$pdf2[x$teta.ind2]), mm(p2pdf2$p2[x$teta.ind2]-p2pdf2$pdf2[x$teta.ind2]), nC2, theta[x$teta.ind2], dof)

  p12[x$teta.ind2] <- C11[x$teta.ind2] - C01[x$teta.ind2] - C10[x$teta.ind2] + C00[x$teta.ind2]
  p12[x$teta.ind2] <- ifelse(p12[x$teta.ind2] < epsilon, epsilon, p12[x$teta.ind2])

                     }


if(length(theta) == 1){  

  C11[x$teta.ind2] <- BiCDF(p1pdf1$p2[x$teta.ind2],                              p2pdf2$p2[x$teta.ind2],                              nC2, theta, dof)
  C01[x$teta.ind2] <- BiCDF(mm(p1pdf1$p2[x$teta.ind2]-p1pdf1$pdf2[x$teta.ind2]), p2pdf2$p2[x$teta.ind2],                              nC2, theta, dof)
  C10[x$teta.ind2] <- BiCDF(p1pdf1$p2[x$teta.ind2],                              mm(p2pdf2$p2[x$teta.ind2]-p2pdf2$pdf2[x$teta.ind2]), nC2, theta, dof)
  C00[x$teta.ind2] <- BiCDF(mm(p1pdf1$p2[x$teta.ind2]-p1pdf1$pdf2[x$teta.ind2]), mm(p2pdf2$p2[x$teta.ind2]-p2pdf2$pdf2[x$teta.ind2]), nC2, theta, dof)

  p12[x$teta.ind2] <- C11[x$teta.ind2] - C01[x$teta.ind2] - C10[x$teta.ind2] + C00[x$teta.ind2]
  p12[x$teta.ind2] <- ifelse(p12[x$teta.ind2] < epsilon, epsilon, p12[x$teta.ind2])

                     } 
                     
                          } 



                                                                                       }                          
                          
                          
                          
                          

}



if(!(x$BivD %in% x$BivD2)){


if(x$margins[1] %in% c(x$VC$m2, x$VC$m3) && x$margins[2] %in% c(x$VC$m2, x$VC$m3))  p12 <- BiCDF(p1, p2, x$nC, theta, dof)


if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){

C1 <- BiCDF(p1pdf1$p2, p2, x$nC, theta, dof) 
C2 <- BiCDF(mm(p1pdf1$p2 - p1pdf1$pdf2), p2, x$nC, theta, dof)  

p12 <- ifelse(C1 - C2 < epsilon, epsilon, C1 - C2)                                  }



if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){

  C11 <- BiCDF(p1pdf1$p2,                 p2pdf2$p2,                 x$nC, theta, dof)
  C01 <- BiCDF(mm(p1pdf1$p2-p1pdf1$pdf2), p2pdf2$p2,                 x$nC, theta, dof)
  C10 <- BiCDF(p1pdf1$p2,                 mm(p2pdf2$p2-p2pdf2$pdf2), x$nC, theta, dof)
  C00 <- BiCDF(mm(p1pdf1$p2-p1pdf1$pdf2), mm(p2pdf2$p2-p2pdf2$pdf2), x$nC, theta, dof)

 p12 <- C11 - C01 - C10 + C00
 p12 <- ifelse(p12 < epsilon, epsilon, p12)  
                                                                                       } 





}




if(cond == 1) p12 <- p12/p1
if(cond == 2) p12 <- p12/p2




if(intervals == TRUE){



bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  
lf <- length(x$coefficients)


#############  
# etas
#############  

# try with 1 number

if(!missing(newdata)){ X1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       X2 <- x$X2 }                       


eta1s <- eta.tr( X1%*%t(bs[,1:x$X1.d2])                     , x$VC$margins[1]) 
eta2s <- eta.tr( X2%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])

#############  
# thetas
#############  



if(  is.null(x$X3) ){ epds <- bs[, lf]
                      if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3) ) dofs <- bs[, lf - 1] else dofs <- dof
                    } 
  
  
  
  
if( !is.null(x$X3) ){ 
  	
  	

 if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){ 
  
  if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont2par){
  
  
       if(!missing(newdata)){
                            X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix") 
                            X6 <- predict.SemiParBIVProbit(x, eq = 6, newdata = newdata, type = "lpmatrix")               
                            }  
       if( missing(newdata)){X5 <- x$X5; X6 <- x$X6} 
 
       dofs <- X5%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
       epds <- X6%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)])
 
                                                                     }
 
  if((x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par) || (x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par) ){
  
       if(!missing(newdata)){
                             X6 <- predict.SemiParBIVProbit(x, eq = 6, newdata = newdata, type = "lpmatrix")  
                             X7 <- predict.SemiParBIVProbit(x, eq = 7, newdata = newdata, type = "lpmatrix")               
                            }
       
       if( missing(newdata)){X6 <- x$X6; X7 <- x$X7}   
  
       dofs <- X6%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)])
       epds <- X7%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2)])
  
                                            }
  
  if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par){
  
       if(!missing(newdata)){
                             X7 <- predict.SemiParBIVProbit(x, eq = 7, newdata = newdata, type = "lpmatrix")
                             X8 <- predict.SemiParBIVProbit(x, eq = 8, newdata = newdata, type = "lpmatrix")               
                             }
       
       if( missing(newdata)){X7 <- x$X7; X8 <- x$X8}    
  
       dofs <- X7%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2)])
       epds <- X8%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2+x$X8.d2)])
  
     }  
  
  
  
                }else{
                
  dofs <- dof   
  
  if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont1par){
  
  
       if(!missing(newdata)) X3 <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X3 <- x$X3 
 
       epds <- X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)])
 
                                                                     }  
  
   if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont2par){
  
  
       if(!missing(newdata)) X4 <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X4 <- x$X4 
 
       epds <- X4%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)])
 
                                                                     }   
  
  
   if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont3par){
  
  
       if(!missing(newdata)) X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X5 <- x$X5 
 
       epds <- X5%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
 
                                                                     }   
  
  
       
       
  if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont2par){
  
  
       if(!missing(newdata)) X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X5 <- x$X5 
 
       epds <- X5%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
 
                                                                     }
 
  if((x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par) || (x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par) ){
  
       if(!missing(newdata)) X6 <- predict.SemiParBIVProbit(x, eq = 6, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X6 <- x$X6   
  
       epds <- X6%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)])
  
                                                                                                                                            }
  
  if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par){
  
       if(!missing(newdata)) X7 <- predict.SemiParBIVProbit(x, eq = 7, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X7 <- x$X7    
  
       epds <- X7%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2)])
  
                                                                    }                
                
                
                
                
                
}  # is.null              
                
                
                
                
                
                
                
                
                


}
  	                        


est.RHOb <- teta.tr(x$VC, epds)$teta
if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3) ) dofs <- dof.tr(dofs)$vao else dofs <- dof 

   
######   
   
if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)) tc <- 1 else tc <- 0  
      
#############  
# sigmas
#############  

      if( is.null(x$X3) ) { # could remove the conditions below and make it more general but need to define
                            # XX.d2 and change the multiplication operator. For now we will leave it like this   
  
if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% c(cont2par,cont3par) ){ ps2 <- lf - (1+tc) }  
if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont2par ){ ps1 <- lf - (2+tc); ps2 <- lf - (1+tc) }
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){ ps1 <- lf - (4+tc); ps2 <- lf - (3+tc) }
if((x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par) || (x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par) ){ ps1 <- lf - (3+tc); ps2 <- lf - (2+tc)}
      
        if(!(x$VC$margins[1] %in% cont1par) ) sigma2.1.star <- bs[, ps1] 
        if(!(x$VC$margins[2] %in% cont1par) ) sigma2.2.star <- bs[, ps2] 
                                
                                } ### is.null
  
  
      if( !is.null(x$X3) ) {


if(x$VC$margins[1] %in% c(cont2par,cont3par) && x$VC$margins[2] %in% c(cont2par,cont3par) )  {    
      
       if(!missing(newdata)){ X3 <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")
                              X4 <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix") }
       if( missing(newdata)){ X3 <- x$X3; X4 <- x$X4 }        
      
       sigma2.1.star <- X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]) 
       sigma2.2.star <- X4%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)]) 

                                                                                              }


if(x$VC$margins[1] %in% c(cont1par) && x$VC$margins[2] %in% c(cont2par,cont3par) ){    
      
       if(!missing(newdata)){ X3 <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix") }
       if( missing(newdata)){ X3 <- x$X3}        
      
       sigma2.2.star <- X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]) 

                                                                                   }

}  ### is.null


    if(!(x$VC$margins[1] %in% cont1par) ) sigma21 <- esp.tr(sigma2.1.star, x$VC$margins[1])$vrb else sigma21 <- 1  
    if(!(x$VC$margins[2] %in% cont1par) ) sigma22 <- esp.tr(sigma2.2.star, x$VC$margins[2])$vrb else sigma22 <- 1   
 
 

 
#############  
# NUs
#############    
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3) )  {    pn1 <- lf - (2+tc) 
                            pn2 <- lf - (1+tc)
                            nu1.st <- bs[, pn1]    # t(as.matrix(bs[, pn1]))
                            nu2.st <- bs[, pn2]  } # t(as.matrix(bs[, pn2]))  } 
  
  if( !is.null(x$X3) ) {  
  
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")
                              X6 <- predict.SemiParBIVProbit(x, eq = 6, newdata = newdata, type = "lpmatrix") }
       if( missing(newdata)){ X5 <- x$X5; X6 <- x$X6 }   
 
       nu1.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
       nu2.st <- X6%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])
  
                       }   
   
nu1 <- esp.tr(nu1.st, x$VC$margins[1])$vrb   
nu2 <- esp.tr(nu2.st, x$VC$margins[2])$vrb   
  
} 


if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  {  pn2 <- lf - (1+tc); nu2.st <- bs[, pn2]  } # nu2.st <- t(as.matrix(bs[, pn2]))
  
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
       if( missing(newdata)){ X5 <- x$X5}     
  
  if( !is.null(x$X3) ) nu2.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
  
  nu2 <- esp.tr(nu2.st, x$VC$margins[2])$vrb   

} 



if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  {  pn2 <- lf - (1+tc); nu2.st <- bs[, pn2]  } # nu2.st <- t(as.matrix(bs[, pn2]))
  
  
       if(!missing(newdata)){ X4 <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")}
       if( missing(newdata)){ X4 <- x$X4}     
  
  if( !is.null(x$X3) ) nu2.st <- X4%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
  
  nu2 <- esp.tr(nu2.st, x$VC$margins[2])$vrb   

}




  
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par ){  
    
  if( is.null(x$X3) )  {  pn1 <- lf - (1+tc); nu1.st <- bs[, pn1]  } # nu1.st <- t(as.matrix(bs[, pn1])) 
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
       if( missing(newdata)){ X5 <- x$X5}    
  
  if( !is.null(x$X3) ) nu1.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
  
  nu1 <- esp.tr(nu1.st, x$VC$margins[1])$vrb  
   
}  


####################################################


if( is.null(x$X3) ){

est.RHOb <- matrix(rep(est.RHOb, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)) dofs <- matrix(rep(dofs, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE) else dofs <- dof

if(is.null(sigma21)) sigma21 <- 1
if(is.null(sigma22)) sigma22 <- 1
if(is.null(nu1))     nu1 <- 1
if(is.null(nu2))     nu2 <- 1

sigma21  <- matrix(rep(sigma21,  each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
sigma22  <- matrix(rep(sigma22,  each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
nu1      <- matrix(rep(nu1,      each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
nu2      <- matrix(rep(nu2,      each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)

                   }



if(x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){

p1s  <- distrHsAT(y1, eta1s, sigma21, nu1, x$margins[1])$p2
p2s  <- distrHsAT(y2, eta2s, sigma22, nu2, x$margins[2])$p2

}


if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){

ppdf1 <- distrHsATDiscr(y1, eta1s, sigma21, nu = 1, x$margins[1], x$VC$y1m) 
p1s   <- ppdf1$p2 
pdf1s <- ppdf1$pdf2

p2s  <- distrHsAT(y2, eta2s, sigma22, nu2, x$margins[2])$p2

}


if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d,x$VC$m2d)){

ppdf1 <- distrHsATDiscr(y1, eta1s, sigma21, nu = 1, x$margins[1], x$VC$y1m) 
p1s   <- ppdf1$p2 
pdf1s <- ppdf1$pdf2

ppdf2 <- distrHsATDiscr(y2, eta2s, sigma22, nu = 1, x$margins[2], x$VC$y2m) 
p2s   <- ppdf2$p2 
pdf2s <- ppdf2$pdf2


}




if(x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){


if(x$VC$BivD %in% c("N","T")) p12s <- matrix(BiCDF(p1s, p2s, x$nC, est.RHOb, dofs, test = FALSE), dim(p1s)[1], n.sim) else{


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(eta1s)[1])

if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- BiCDF(p1s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,])                  
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- BiCDF(p1s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,])
                      
                        }

if(!(x$BivD %in% x$BivD2)) p12s <- BiCDF(p1s, p2s, x$nC, est.RHOb, dofs, test = FALSE)
                                                                                                                          }


                                                                                 }


if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){


if(x$VC$BivD %in% c("N","T")){ C1s <- matrix(BiCDF(p1s, p2s,       x$nC, est.RHOb, dofs, test = FALSE), dim(p1s)[1], n.sim)
                               C2s <- matrix(BiCDF(p1s-pdf1s, p2s, x$nC, est.RHOb, dofs, test = FALSE), dim(p1s)[1], n.sim)
                               p12s <- ifelse(C1s - C2s < epsilon, epsilon, C1s - C2s)

                             }else{


if(x$BivD %in% x$BivD2){

p12s <- C1s <- C2s <- matrix(NA, ncol = n.sim, nrow = dim(eta1s)[1])

if( length(x$teta1) != 0){

C1s[x$teta.ind1,]  <- BiCDF(p1s[x$teta.ind1,],                     p2s[x$teta.ind1,], nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE)
C2s[x$teta.ind1,]  <- BiCDF(p1s[x$teta.ind1,]-pdf1s[x$teta.ind1,], p2s[x$teta.ind1,], nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE)
p12s[x$teta.ind1,] <- ifelse(C1s[x$teta.ind1,] - C2s[x$teta.ind1,] < epsilon, epsilon, C1s[x$teta.ind1,] - C2s[x$teta.ind1,])

                         }                            
 
if( length(x$teta2) != 0){

C1s[x$teta.ind2,]  <- BiCDF(p1s[x$teta.ind2,],                     p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE)
C2s[x$teta.ind2,]  <- BiCDF(p1s[x$teta.ind2,]-pdf1s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE)
p12s[x$teta.ind2,] <- ifelse(C1s[x$teta.ind2,] - C2s[x$teta.ind2,] < epsilon, epsilon, C1s[x$teta.ind2,] - C2s[x$teta.ind2,])

                         }  
                            
                            
                      
                        }



if(!(x$BivD %in% x$BivD2)){ C1s <- BiCDF(p1s,       p2s, x$nC, est.RHOb, dofs, test = FALSE)
                            C2s <- BiCDF(p1s-pdf1s, p2s, x$nC, est.RHOb, dofs, test = FALSE)
                            p12s <- ifelse(C1s - C2s < epsilon, epsilon, C1s - C2s)
                          }


                                   }


                                                                                  }











if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d,x$VC$m2d)){


if(x$VC$BivD %in% c("N","T")){ 

  C11s <- matrix(BiCDF(p1s,       p2s,       x$nC, est.RHOb, dofs, test = FALSE), dim(p1s)[1], n.sim)
  C01s <- matrix(BiCDF(p1s-pdf1s, p2s,       x$nC, est.RHOb, dofs, test = FALSE), dim(p1s)[1], n.sim)
  C10s <- matrix(BiCDF(p1s, p2s-pdf2s,       x$nC, est.RHOb, dofs, test = FALSE), dim(p1s)[1], n.sim)
  C00s <- matrix(BiCDF(p1s-pdf1s, p2s-pdf2s, x$nC, est.RHOb, dofs, test = FALSE), dim(p1s)[1], n.sim)

 p12s <- C11s - C01s - C10s + C00s

                              }else{


if(x$BivD %in% x$BivD2){

p12s <- C11s <- C01s <- C10s <- C00s <- matrix(NA, ncol = n.sim, nrow = dim(eta1s)[1])

if( length(x$teta1) != 0){

  C11s[x$teta.ind1,] <- BiCDF(p1s[x$teta.ind1,],                     p2s[x$teta.ind1,],                     nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE)
  C01s[x$teta.ind1,] <- BiCDF(p1s[x$teta.ind1,]-pdf1s[x$teta.ind1,], p2s[x$teta.ind1,],                     nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE)
  C10s[x$teta.ind1,] <- BiCDF(p1s[x$teta.ind1,],                     p2s[x$teta.ind1,]-pdf2s[x$teta.ind1,], nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE)
  C00s[x$teta.ind1,] <- BiCDF(p1s[x$teta.ind1,]-pdf1s[x$teta.ind1,], p2s[x$teta.ind1,]-pdf2s[x$teta.ind1,], nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE)

 p12s[x$teta.ind1,] <- C11s[x$teta.ind1,] - C01s[x$teta.ind1,] - C10s[x$teta.ind1,] + C00s[x$teta.ind1,]  
 p12s[x$teta.ind1,] <- ifelse(p12s[x$teta.ind1,] < epsilon, epsilon, p12s[x$teta.ind1,]) 

}


if( length(x$teta2) != 0){

  C11s[x$teta.ind2,] <- BiCDF(p1s[x$teta.ind2,],                     p2s[x$teta.ind2,],                     nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE)
  C01s[x$teta.ind2,] <- BiCDF(p1s[x$teta.ind2,]-pdf1s[x$teta.ind2,], p2s[x$teta.ind2,],                     nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE)
  C10s[x$teta.ind2,] <- BiCDF(p1s[x$teta.ind2,],                     p2s[x$teta.ind2,]-pdf2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE)
  C00s[x$teta.ind2,] <- BiCDF(p1s[x$teta.ind2,]-pdf1s[x$teta.ind2,], p2s[x$teta.ind2,]-pdf2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE)

 p12s[x$teta.ind2,] <- C11s[x$teta.ind2,] - C01s[x$teta.ind2,] - C10s[x$teta.ind2,] + C00s[x$teta.ind2,]  
 p12s[x$teta.ind2,] <- ifelse(p12s[x$teta.ind2,] < epsilon, epsilon, p12s[x$teta.ind2,]) 

}



                        }



if(!(x$BivD %in% x$BivD2)){

                           
  C11s <- BiCDF(p1s, p2s,             x$nC, est.RHOb, dofs, test = FALSE)
  C01s <- BiCDF(p1s-pdf1s, p2s,       x$nC, est.RHOb, dofs, test = FALSE)
  C10s <- BiCDF(p1s, p2s-pdf2s,       x$nC, est.RHOb, dofs, test = FALSE)
  C00s <- BiCDF(p1s-pdf1s, p2s-pdf2s, x$nC, est.RHOb, dofs, test = FALSE)

 p12s <- C11s - C01s - C10s + C00s  
 p12s <- ifelse(p12s < epsilon, epsilon, p12s)  

                           
                           


                          }


                                     }


}














if(cond == 1) p12s <- p12s/p1s
if(cond == 2) p12s <- p12s/p2s





} # interv






                        }## biv
                        
######################################################################################################
######################################################################################################





######################################################################################################
######################################################################################################

      if(type == "independence"){



if(!missing(newdata)){

nu1 <- nu2 <- sigma21 <- sigma22 <- NA

eta1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$coefficients[1:x$X1.d2]
eta2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$coefficients[1:x$X2.d2]

if( !is.null(x$X3) ){ 


if(x$margins[1] %in% cont1par && x$margins[2] %in% c(cont2par,cont3par))             sigma22 <- esp.tr(predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$coefficients[(x$X2.d2+1):(x$X2.d2+x$X3.d2)], x$margins[2])$vrb
if(x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par,cont3par)) sigma22 <- esp.tr(predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$coefficients[(x$X2.d2+1):(x$X2.d2+x$X4.d2)], x$margins[2])$vrb
if(!(x$margins[1] %in% cont1par))                                                    sigma21 <- esp.tr(predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$coefficients[(x$X1.d2+1):(x$X1.d2+x$X3.d2)], x$margins[1])$vrb


if(x$margins[1] %in% cont1par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 4  ; indn1 <- NA; indn2 <- (x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){ eq.nu1 <- 5;  eq.nu2 <- 6  ; indn1 <- (x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2); indn2 <- (x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X6.d2)}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA ; indn1 <- NA; indn2 <- NA}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 5  ; indn1 <- NA; indn2 <- (x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X5.d2)}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){ eq.nu1 <- 5;  eq.nu2 <- NA ; indn1 <- (x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2); indn2 <- NA}

if(x$margins[1] %in% cont3par) nu1 <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu1, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$coefficients[indn1], x$margins[1])$vrb
if(x$margins[2] %in% cont3par) nu2 <- esp.tr(predict.SemiParBIVProbit(x, eq = eq.nu2, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$coefficients[indn2], x$margins[2])$vrb

}

if( is.null(x$X3) ){ 

sigma21 <- x$gamlss1$sigma2
sigma22 <- x$gamlss2$sigma2

nu1 <- x$gamlss1$nu
nu2 <- x$gamlss2$nu

}

}



if(missing(newdata)){


eta1 <- x$gamlss1$eta1
eta2 <- x$gamlss2$eta1

sigma21 <- x$gamlss1$sigma2
sigma22 <- x$gamlss2$sigma2

nu1 <- x$gamlss1$nu
nu2 <- x$gamlss2$nu

}






if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p1 <- distrHsAT(y1, eta1, sigma21, nu1, x$margins[1])$p2 
if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p2 <- distrHsAT(y2, eta2, sigma22, nu2, x$margins[2])$p2 

if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d)){ 
                                           p1pdf1 <- distrHsATDiscr(y1, eta1, sigma21, nu = 1, x$margins[1], x$VC$y1m)
                                           p1     <- p1pdf1$p2 - p1pdf1$pdf2
                                           }
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){ 
                                           p2pdf2 <- distrHsATDiscr(y2, eta2, sigma22, nu = 1, x$margins[2], x$VC$y2m)
                                           p2     <- p2pdf2$p2 - p2pdf2$pdf2
                                           }  

p12 <- p1*p2

if(cond == 1) p12 <- p2
if(cond == 2) p12 <- p1



if(intervals == TRUE){

bs1 <- rMVN(n.sim, mean = x$gamlss1$coefficients, sigma=x$gamlss1$Vb)
bs2 <- rMVN(n.sim, mean = x$gamlss2$coefficients, sigma=x$gamlss2$Vb)


#############  
# etas
############# 

if(!missing(newdata)){ X1 <- predict.SemiParBIVProbit(x, eq = 1, newdata = newdata, type = "lpmatrix")
                       X2 <- predict.SemiParBIVProbit(x, eq = 2, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X1 <- x$X1; X2 <- x$X2}  

eta1s <- eta.tr( X1%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1]) 
eta2s <- eta.tr( X2%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2]) 

#############  
# sigmas
#############  

      if( is.null(x$X3) ) {
      
      
      if(x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par,cont3par)){
      
       sigma2.1.star <- bs1[, x$X1.d2 + 1] 
       sigma2.2.star <- bs2[, x$X2.d2 + 1] 
      
      }
      
      if(x$margins[1] %in% c(cont1par) && x$margins[2] %in% c(cont2par,cont3par)){
       
       sigma2.2.star <- bs2[, x$X2.d2 + 1] 
      
      }      
      
      
                          }
  
  
  
  
      if( !is.null(x$X3) ) {

if(x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par,cont3par)){

if(!missing(newdata)){ X3 <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")
                       X4 <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X3 <- x$X3; X4 <- x$X4}  

       sigma2.1.star <- X3%*%t(bs1[,(x$X1.d2+1):(x$X1.d2+x$X3.d2)]) 
       sigma2.2.star <- X4%*%t(bs2[,(x$X2.d2+1):(x$X2.d2+x$X4.d2)]) 
                                                                                     }
                                                                                     
if(x$margins[1] %in% c(cont1par) && x$margins[2] %in% c(cont2par,cont3par)){

if(!missing(newdata)) X3 <- predict.SemiParBIVProbit(x, eq = 3, newdata = newdata, type = "lpmatrix")  
if( missing(newdata)) X3 <- x$X3 

       sigma2.2.star <- X3%*%t(bs2[,(x$X2.d2+1):(x$X2.d2+x$X3.d2)]) 
                                                                            }                                                                                     
                           }  

if(!(x$margins[1] %in% c(cont1par))) sigma21 <- esp.tr(sigma2.1.star, x$VC$margins[1])$vrb   
                                     sigma22 <- esp.tr(sigma2.2.star, x$VC$margins[2])$vrb  

#############  
# NUs
#############    
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3) )  {     
       nu1.st <- bs1[, x$X1.d2 + 2]  # t(as.matrix(bs1[, x$X1.d2 + 2]))
       nu2.st <- bs2[, x$X2.d2 + 2]  # t(as.matrix(bs2[, x$X2.d2 + 2]))   
                        } 
  
  if( !is.null(x$X3) ) {  
  
  
if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")
                       X6 <- predict.SemiParBIVProbit(x, eq = 6, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X5 <- x$X5; X6 <- x$X6}    
  
       nu1.st <- X5%*%t(bs1[,(x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2)]) 
       nu2.st <- X6%*%t(bs2[,(x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X6.d2)])
                       }   
   
nu1 <- esp.tr(nu1.st, x$VC$margins[1])$vrb   
nu2 <- esp.tr(nu2.st, x$VC$margins[2])$vrb   
  
} 


if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  nu2.st <- bs2[, x$X2.d2 + 2] # t(as.matrix(bs2[, x$X2.d2 + 2]))  
  
  
  if( !is.null(x$X3) ){
  
  
    if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
    if( missing(newdata)){ X5 <- x$X5}  
  
    nu2.st <- X5%*%t(bs2[,(x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X5.d2)]) 
  
  
                      }
  
  
                       nu2    <- esp.tr(nu2.st, x$VC$margins[2])$vrb   
} 



if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  nu2.st <- bs2[, x$X2.d2 + 1]   
  
  
  if( !is.null(x$X3) ){
  
  
    if(!missing(newdata)){ X4 <- predict.SemiParBIVProbit(x, eq = 4, newdata = newdata, type = "lpmatrix")}
    if( missing(newdata)){ X4 <- x$X4}  
  
    nu2.st <- X4%*%t(bs2[,(x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
  
  
                      }
  
   nu2    <- esp.tr(nu2.st, x$VC$margins[2])$vrb   
} 


  
  
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par ){  
    
  if( is.null(x$X3) )  nu1.st <- bs1[, x$X1.d2 + 2] # t(as.matrix(bs1[, x$X1.d2 + 2]))   


  if( !is.null(x$X3) ){
    
    if(!missing(newdata)){ X5 <- predict.SemiParBIVProbit(x, eq = 5, newdata = newdata, type = "lpmatrix")}
    if( missing(newdata)){ X5 <- x$X5}    
  
    nu1.st <- X5%*%t(bs1[,(x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2)])                  
  
                      }
  
  nu1    <- esp.tr(nu1.st, x$VC$margins[1])$vrb    
}  


if( is.null(x$X3) ){

sigma21  <- matrix(rep(sigma21, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
sigma22  <- matrix(rep(sigma22, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
nu1      <- matrix(rep(nu1, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
nu2      <- matrix(rep(nu2, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)

                   }




if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p1s <- distrHsAT(y1, eta1s, sigma21, nu1, x$margins[1])$p2 
if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p2s <- distrHsAT(y2, eta2s, sigma22, nu2, x$margins[2])$p2

if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d)){ 
                                           p1pdf1s <- distrHsATDiscr(y1, eta1s, sigma21, nu = 1, x$margins[1], x$VC$y1m)
                                           p1s     <- p1pdf1s$p2 - p1pdf1s$pdf2
                                           }
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){ 
                                           p2pdf2s <- distrHsATDiscr(y2, eta2s, sigma22, nu = 1, x$margins[2], x$VC$y2m)
                                           p2s     <- p2pdf2s$p2 - p2pdf2s$pdf2
                                           }  

p12s <- p1s*p2s

if(cond == 1) p12s <- p2s
if(cond == 2) p12s <- p1s



} # intervals


} # independence



list(p12 = p12, p12s = p12s, p1 = p1, p2 = p2)


}




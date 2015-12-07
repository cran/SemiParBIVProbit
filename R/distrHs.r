distrHs <- function(y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2, naive = FALSE){


p2 <- derp2.dersigma.st <- derp2.dereta2 <- der2p2.dereta2eta2 <- der2p2.dersigma2.st2 <- der2p2.dereta2dersigma2.st <- 1

der2pdf2.dereta2dernu.st    = 1
der2pdf2.sigma2.st2dernu.st = 1
derpdf2.dernu.st            = 1
der2pdf2.dernu.st2          = 1
derp2.nu.st                 = 1
der2p2.dernu.st2            = 1
der2p2.dereta2dernu.st      = 1
der2p2.dersigma2.stdernu.st = 1



eps <- 1e-06; eps2 <- eps*eps 

cont2par <- c("WEI","WEI2","iG","LO","rGU","GU","GA","iGA")
cont3par <- c("DAGUM")





if(naive == FALSE) {

if(margin2 %in% c("GA"))         fp <- function(sigma2) pgamma(y2, shape = 1/sigma2, scale = exp(eta2) * sigma2)
if(margin2 == "iGA")             fp <- function(sigma2) 1 - pgamma(((exp(eta2) * (1/sigma2 + 1))/y2), shape = 1/sigma2, scale=1)


if(margin2 %in% c("GA","iGA")){

sigma21 <- sigma2 - eps/2
sigma22 <- sigma2 + eps/2 
f1 <- fp(sigma21)
f2 <- fp(sigma22)

derp2.dersigma2 <- (f2 - f1)/ eps   

    
t01 <- t10 <- t11 <- t01 <- sigma2 + eps
t11 <- t11 + eps 

f00 <- fp(sigma2) 
f01 <- fp(t01) 
f10 <- fp(t10) 
f11 <- fp(t11) 

der2p2.dersigma22 <- (f11 - f01 - f10 + f00)/eps2

     }

}





if(margin2 %in% c("N","LN")){

  pdf2          <- dnorm(y2, mean=eta2, sd = sqrt(sigma2))
derpdf2.dereta2 <- (1/(sqrt(2 * pi * sigma2))) * (exp(-0.5 * (y2 - eta2)^2/sigma2) * (0.5 * (2 * (y2 - eta2))/sigma2))     
derpdf2.sigma2 <-  (-((1/(sqrt(2 * pi * sigma2))) * (exp(-0.5 * (y2 - eta2)^2/sigma2) * 
                   (-0.5 * (y2 - eta2)^2/sigma2^2)) + 0.5 * (2 * pi * (2 * pi * 
                   sigma2)^-0.5)/(sqrt(2 * pi * sigma2))^2 * exp(-0.5 * (y2 - 
                   eta2)^2/sigma2)))  
dersigma2.dersigma2.st <- exp(sigma2.st)   
derpdf2.dersigma2.st <- derpdf2.sigma2 * dersigma2.dersigma2.st   # left here for a reason
der2pdf2.dereta2 <-  (1/(sqrt(2 * pi * sigma2))) * (exp(-0.5 * (y2 - eta2)^2/sigma2) * 
                     (0.5 * (2 * (y2 - eta2))/sigma2) * (0.5 * (2 * (y2 - eta2))/sigma2) - 
                     exp(-0.5 * (y2 - eta2)^2/sigma2) * (0.5 * 2/sigma2))  
                     
                     
der2pdf2.dersigma2.st2 <- (0.5 * (2 * pi * (2 * pi * exp(sigma2.st))^-0.5)/(sqrt(2 * pi * 
                          exp(sigma2.st)))^2 * (exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) * 
                          (-0.5 * (y2 - eta2)^2 * exp(sigma2.st)/exp(sigma2.st)^2)) + 
                          (0.5 * (2 * pi * ((2 * pi * exp(sigma2.st))^-(0.5 + 1) * 
                          (0.5 * (2 * pi * exp(sigma2.st)))))/(sqrt(2 * pi * exp(sigma2.st)))^2 + 
                          0.5 * (2 * pi * (2 * pi * exp(sigma2.st))^-0.5) * (2 * 
                          (0.5 * (2 * pi * exp(sigma2.st) * (2 * pi * exp(sigma2.st))^-0.5) * 
                          (sqrt(2 * pi * exp(sigma2.st)))))/((sqrt(2 * 
                          pi * exp(sigma2.st)))^2)^2) * exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) + 
                          ((1/(sqrt(2 * pi * exp(sigma2.st)))) * (exp(-0.5 * (y2 - 
                          eta2)^2/exp(sigma2.st)) * (-0.5 * (y2 - eta2)^2 * (2 * 
                          (exp(sigma2.st) * exp(sigma2.st)))/(exp(sigma2.st)^2)^2) + 
                          exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) * (-0.5 * (y2 - 
                          eta2)^2 * exp(sigma2.st)/exp(sigma2.st)^2) * (-0.5 * 
                          (y2 - eta2)^2/exp(sigma2.st)^2)) + 0.5 * (2 * pi * 
                          exp(sigma2.st) * (2 * pi * exp(sigma2.st))^-0.5)/(sqrt(2 * 
                          pi * exp(sigma2.st)))^2 * (exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) * 
                          (-0.5 * (y2 - eta2)^2/exp(sigma2.st)^2)))) * exp(sigma2.st) + 
                          (-((1/(sqrt(2 * pi * exp(sigma2.st)))) * (exp(-0.5 * (y2 - 
                          eta2)^2/exp(sigma2.st)) * (-0.5 * (y2 - eta2)^2/exp(sigma2.st)^2)) + 
                          0.5 * (2 * pi * (2 * pi * exp(sigma2.st))^-0.5)/(sqrt(2 * 
                          pi * exp(sigma2.st)))^2 * exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)))) * 
                          exp(sigma2.st)  
                          
   
   
der2pdf2.dereta2dersigma2.st <- -((1/(sqrt(2 * pi * exp(sigma2.st)))) * (exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) * 
                                 (0.5 * (2 * (y2 - eta2)) * exp(sigma2.st)/exp(sigma2.st)^2) + 
                                 exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) * (-0.5 * (y2 - 
                                 eta2)^2 * exp(sigma2.st)/exp(sigma2.st)^2) * (0.5 * (2 * 
                                 (y2 - eta2))/exp(sigma2.st))) + 0.5 * (2 * pi * exp(sigma2.st) * 
                                 (2 * pi * exp(sigma2.st))^-0.5)/(sqrt(2 * pi * exp(sigma2.st)))^2 * 
                                 (exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) * (0.5 * (2 * (y2 - 
                                 eta2))/exp(sigma2.st))))  
  
  
  
if(naive == FALSE){  
  
    p2          <- pnorm(y2, mean=eta2, sd = sqrt(sigma2))

derp2.dereta2    <- -pdf2
                    
derp2.dersigma.st <- 0.5*(y2-eta2)*derp2.dereta2     


der2p2.dereta2eta2 <- -((1/(sqrt(2 * pi * sigma2))) * (exp(-0.5 * (y2 - eta2)^2/sigma2) * 
                        (0.5 * (2 * (y2 - eta2))/sigma2)))


der2p2.dersigma2.st2 <-  0.5 * (y2 - eta2) * (0.5 * (2 * pi * exp(sigma2.st) * (2 * pi * 
                         exp(sigma2.st))^-0.5)/(sqrt(2 * pi * exp(sigma2.st)))^2 * 
                         exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) + (1/(sqrt(2 * pi * 
                         exp(sigma2.st)))) * (exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) * 
                         (-0.5 * (y2 - eta2)^2 * exp(sigma2.st)/exp(sigma2.st)^2)))

der2p2.dereta2dersigma2.st <-  0.5 * (2 * pi * exp(sigma2.st) * (2 * pi * exp(sigma2.st))^-0.5)/(sqrt(2 * 
                               pi * exp(sigma2.st)))^2 * exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) + 
                               (1/(sqrt(2 * pi * exp(sigma2.st)))) * (exp(-0.5 * (y2 - eta2)^2/exp(sigma2.st)) * 
                               (-0.5 * (y2 - eta2)^2 * exp(sigma2.st)/exp(sigma2.st)^2))          
}


}

####


if(margin2 == "DAGUM"){




#a <- sqrt(sigma2)  shape 1 
#p <- nu            shape 2 
#b <- exp(eta2)     scale                          



pdf2 <- sqrt(sigma2)*nu/y2*( ((y2/exp(eta2))^(sqrt(sigma2)*nu))/  ( (y2/exp(eta2))^sqrt(sigma2) + 1 )^(nu+1) )            
  
  
derpdf2.dereta2 <- (y2^(nu*sqrt(sigma2))*(nu*sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^2*sigma2*exp(2*eta2*sqrt(sigma2))))/((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*(y2^(2*sqrt(sigma2)+1)+2*exp(eta2*sqrt(sigma2))*y2^(sqrt(sigma2)+1)+exp(2*eta2*sqrt(sigma2))*y2))    
    
                    
derpdf2.sigma2 <- -(sqrt(sigma2)*y2^(nu*sqrt(sigma2))*(nu*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^2*exp(2*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2)+y2^(nu*sqrt(sigma2))*(-nu*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu*exp(2*eta2*sqrt(sigma2))))/(sqrt(sigma2)*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*(2*y2^(2*sqrt(sigma2)+1)+4*exp(eta2*sqrt(sigma2))*y2^(sqrt(sigma2)+1)+2*exp(2*eta2*sqrt(sigma2))*y2))    
    

derpdf2.nu <- -(sqrt(sigma2)*(nu*exp(eta2*sqrt(sigma2))*y2^(nu*sqrt(sigma2))*log(exp(-eta2*sqrt(sigma2))*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2))))-exp(eta2*sqrt(sigma2))*y2^(nu*sqrt(sigma2)))-nu*sigma2*
exp(eta2*sqrt(sigma2))*y2^(nu*sqrt(sigma2))*log(exp(-eta2)*y2))/((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*(y2^(sqrt(sigma2)+1)+exp(eta2*sqrt(sigma2))*y2))


dersigma2.dersigma2.st <- exp(sigma2.st)   

dernu.dernu.st <- exp(nu.st)


 der2pdf2.dereta2 <- (sqrt(sigma2)*y2^(nu*sqrt(sigma2))*(nu*sigma2*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+(-3*nu^2-nu)*sigma2*exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2)+nu^3*sigma2*exp(3*eta2*sqrt(sigma2))))/((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*(y2^(3*sqrt(sigma2)+1)+3*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2)+1)+3*exp(2*eta2*sqrt(sigma2))*y2^(sqrt(sigma2)+1)+exp(3*eta2*sqrt(sigma2))*y2))  
  
  
der2pdf2.dersigma22 <- -(y2^(nu*sqrt(sigma2))*(-nu*sigma2*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+(3*nu^2+nu)*sigma2*exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^3*sigma2*exp(3*eta2*sqrt(sigma2)))*
log(exp(-eta2)*y2)^2+sqrt(sigma2)*y2^(nu*sqrt(sigma2))*(nu*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+(nu-nu^2)*exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^2*exp(3*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2)+
y2^(nu*sqrt(sigma2))*(nu*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+2*nu*exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2)+nu*exp(3*eta2*sqrt(sigma2))))/(sqrt(sigma2)*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*
(4*sigma2*y2^(3*sqrt(sigma2)+1)+12*sigma2*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2)+1)+12*sigma2*exp(2*eta2*sqrt(sigma2))*y2^(sqrt(sigma2)+1)+4*sigma2*exp(3*eta2*sqrt(sigma2))*y2))




    
    der2pdf2.dernu2 <- (sqrt(sigma2)*(nu*exp(eta2*sqrt(sigma2))*y2^(nu*sqrt(sigma2))*log(exp(-eta2*sqrt(sigma2))*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2))))^2-2*exp(eta2*sqrt(sigma2))*y2^(nu*sqrt(sigma2))*
log(exp(-eta2*sqrt(sigma2))*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2))))+nu*sigma2*exp(eta2*sqrt(sigma2))*y2^(nu*sqrt(sigma2))*log(exp(-eta2)*y2)^2)-2*nu*sigma2*exp(eta2*sqrt(sigma2))*y2^(nu*sqrt(sigma2))*
log(exp(-eta2)*y2)*log(exp(-eta2*sqrt(sigma2))*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2))))+2*sigma2*exp(eta2*sqrt(sigma2))*y2^(nu*sqrt(sigma2))*log(exp(-eta2)*y2))/(
(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*(y2^(sqrt(sigma2)+1)+exp(eta2*sqrt(sigma2))*y2))   
    
    
    
    
der2pdf2.dereta2dersigma2 <- -(sqrt(sigma2)*y2^(nu*sqrt(sigma2))*(nu*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+(-3*nu^2-nu)*exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2)+nu^3*exp(3*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2)+y2^(nu*sqrt(sigma2))*
 (-2*nu*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+(2*nu^2-2*nu)*exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2)+2*nu^2*exp(3*eta2*sqrt(sigma2))))/((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*
(2*y2^(3*sqrt(sigma2)+1)+6*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2)+1)+6*exp(2*eta2*sqrt(sigma2))*y2^(sqrt(sigma2)+1)+2*exp(3*eta2*sqrt(sigma2))*y2))
 
    
 der2pdf2.dereta2dernu <- (y2^(nu*sqrt(sigma2))*(nu^2*sigma2*exp(2*eta2*sqrt(sigma2))-nu*sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*log(exp(-eta2*sqrt(sigma2))*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2))))+sqrt(sigma2)*
y2^(nu*sqrt(sigma2))*(nu*sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^2*sigma2*exp(2*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2)+y2^(nu*sqrt(sigma2))*
(sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-2*nu*sigma2*exp(2*eta2*sqrt(sigma2))))/((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*
(y2^(2*sqrt(sigma2)+1)+2*exp(eta2*sqrt(sigma2))*y2^(sqrt(sigma2)+1)+exp(2*eta2*sqrt(sigma2))*y2))    
    
    
    
    der2pdf2.dersigma2dernu <- (sqrt(sigma2)*(y2^(nu*sqrt(sigma2))*(nu*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^2*exp(2*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2)*log(exp(-eta2*sqrt(sigma2))*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2))))+
y2^(nu*sqrt(sigma2))*((nu-1)*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)+3*nu*exp(2*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2))+y2^(nu*sqrt(sigma2))*(-nu*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu*exp(2*eta2*sqrt(sigma2)))*
log(exp(-eta2*sqrt(sigma2))*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2))))+y2^(nu*sqrt(sigma2))*(nu^2*sigma2*exp(2*eta2*sqrt(sigma2))-nu*sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*log(exp(-eta2)*y2)^2+
y2^(nu*sqrt(sigma2))*(exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)+exp(2*eta2*sqrt(sigma2))))/(sqrt(sigma2)*(y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))^nu*
(2*y2^(2*sqrt(sigma2)+1)+4*exp(eta2*sqrt(sigma2))*y2^(sqrt(sigma2)+1)+2*exp(2*eta2*sqrt(sigma2))*y2)) 
 
 
 
 
if(naive == FALSE){   
 
 
 #a <- sqrt(sigma2)
 #p <- nu
 #b <- exp(eta2)         

 
 
    
    p2  <- ( 1 + (y2/exp(eta2))^-sqrt(sigma2) )^-nu 


derp2.dereta2 <- -((1 + (y2/exp(eta2))^-sqrt(sigma2))^-(nu + 1) * (nu * ((y2/exp(eta2))^-(sqrt(sigma2) + 
    1) * (sqrt(sigma2) * (y2 * exp(eta2)/exp(eta2)^2)))))
                
                          
derp2.dersigma2 <- (1 + (y2/exp(eta2))^-sqrt(sigma2))^-(nu + 1) * (nu * ((y2/exp(eta2))^-sqrt(sigma2) * 
    (log((y2/exp(eta2))) * (0.5 * sigma2^-0.5))))           
    

derp2.dernu <- -((1 + (y2/exp(eta2))^-sqrt(sigma2))^-nu * log((1 + (y2/exp(eta2))^-sqrt(sigma2))))


der2p2.dereta2eta2 <- -(nu*sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^2*sigma2*exp(2*eta2*sqrt(sigma2)))/(((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))/y2^sqrt(sigma2))^nu*(y2^(2*sqrt(sigma2))+2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)+exp(2*eta2*sqrt(sigma2))))




der2p2.dersigma22 <- -(sqrt(sigma2)*(nu*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^2*exp(2*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2)^2+(nu*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)+nu*exp(2*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2))/(sqrt(sigma2)*((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))/y2^sqrt(sigma2))^nu*(4*sigma2*y2^(2*sqrt(sigma2))+8*sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)+4*sigma2*exp(2*eta2*sqrt(sigma2))))   
   
   
   
        der2p2.dernu2 <- (1 + (y2/exp(eta2))^-sqrt(sigma2))^-nu * log((1 + (y2/exp(eta2))^-sqrt(sigma2))) * 
    log((1 + (y2/exp(eta2))^-sqrt(sigma2))) 
        
        

der2p2.dereta2dersigma2 <- (sqrt(sigma2)*(nu*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu^2*exp(2*eta2*sqrt(sigma2)))*log(exp(-eta2)*y2)-nu*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)-nu*exp(2*eta2*sqrt(sigma2)))/(sqrt(sigma2)*((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))/y2^sqrt(sigma2))^nu*(2*y2^(2*sqrt(sigma2))+4*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)+2*exp(2*eta2*sqrt(sigma2))))



der2p2.dereta2dernu <- (sqrt(sigma2)*(nu*exp(eta2*sqrt(sigma2))*log((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))/y2^sqrt(sigma2))-exp(eta2*sqrt(sigma2))))/((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))*((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))/y2^sqrt(sigma2))^nu)


der2p2.dersigma2dernu <- -(sqrt(sigma2)*(nu*exp(eta2*sqrt(sigma2))*log(exp(-eta2)*y2)*log((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))/y2^sqrt(sigma2))-exp(eta2*sqrt(sigma2))*log(exp(-eta2)*y2)))/(((y2^sqrt(sigma2)+exp(eta2*sqrt(sigma2)))/y2^sqrt(sigma2))^nu*(2*sigma2*y2^sqrt(sigma2)+2*sigma2*exp(eta2*sqrt(sigma2))))  
  

}






}



####





















####################


if(margin2 == "WEI"){

  pdf2          <- sqrt(sigma2)/exp(eta2)*(y2/exp(eta2))^(sqrt(sigma2)-1) * exp(-(y2/exp(eta2))^sqrt(sigma2))
  
  
derpdf2.dereta2 <- (exp(-(exp(-eta2)* y2)^sqrt(sigma2))* sigma2 *(exp(-eta2)* y2)^sqrt(sigma2)* (-1 + (exp(-eta2)* y2)^sqrt(
                    sigma2)))/y2    
                    
                    
derpdf2.sigma2 <- -(exp(-exp(-eta2*sqrt(sigma2))*y2^sqrt(sigma2)-2*eta2*sqrt(sigma2))*(sqrt(sigma2)*(y2^(2*sqrt(sigma2))-exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*log(exp(-eta2)*y2)-exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)))/(2*sqrt(sigma2)*y2)


dersigma2.dersigma2.st <- exp(sigma2.st)   


 der2pdf2.dereta2 <-  (exp(-(exp(-eta2)* y2)^sqrt(sigma2))* sigma2^(
 3/2)* (exp(-eta2)* y2)^sqrt(sigma2)* (1 - 
   3 *(exp(-eta2)* y2)^sqrt(sigma2) + (exp(-eta2)* y2)^(2* sqrt(sigma2))))/y2    
    
    
der2pdf2.dersigma22 <- -(exp(-exp(-eta2*sqrt(sigma2))*y2^sqrt(sigma2)-3*eta2*sqrt(sigma2))*((-sigma2*y2^(3*sqrt(sigma2))+3*sigma2*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))-sigma2*exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2))*
log(exp(-eta2)*y2)^2+sqrt(sigma2)*(exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))-exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2))*log(exp(-eta2)*y2)+exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2)))/(4*sigma2^(3/2)*y2)    
    
    
    
der2pdf2.dereta2dersigma2 <- -(exp(-exp(-eta2*sqrt(sigma2))*y2^sqrt(sigma2)-3*eta2*sqrt(sigma2))*(sqrt(sigma2)*(y2^(3*sqrt(sigma2))-3*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2))*log(exp(-eta2)*y2)-
2*exp(eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+2*exp(2*eta2*sqrt(sigma2))*y2^sqrt(sigma2)))/(2*y2) 
 
 
 
 
if(naive == FALSE){   
 
    
    p2  <-  1-exp(-(y2/exp(eta2))^sqrt(sigma2)) #  pWEI(y2,exp(eta2), sqrt(sigma2))


derp2.dereta2    <- -(exp(-(y2/exp(eta2))^sqrt(sigma2)) * ((y2/exp(eta2))^(sqrt(sigma2) - 
                      1) * (sqrt(sigma2) * (y2 * exp(eta2)/exp(eta2)^2))))

                
                          
derp2.dersigma2 <-  exp(-(y2/exp(eta2))^sqrt(sigma2)) * ((y2/exp(eta2))^sqrt(sigma2) * 
                   (log((y2/exp(eta2))) * (0.5 * sigma2^-0.5)))               
    


der2p2.dereta2eta2 <- -exp(-(exp(-eta2)* y2)^sqrt(
   sigma2))* sigma2* (exp(-eta2)* y2)^sqrt(sigma2)* (-1 + (exp(-eta2)* y2)^sqrt(
   sigma2))



der2p2.dersigma22 <- (exp(-(exp(-eta2)* y2)^sqrt(sigma2))* (exp(-eta2)* y2)^sqrt(sigma2)*
  log(exp(-eta2)* y2)* (-0.25* sigma2^1 + 
   sigma2^1.5 *(0.25 - 0.25* (exp(-eta2)* y2)^sqrt(sigma2)) *log(
     exp(-eta2)* y2)))/sigma2^2.5


der2p2.dereta2dersigma2 <- (exp(-(exp(-eta2)* y2)^sqrt(
  sigma2))* (exp(-eta2)* y2)^sqrt(sigma2)* (-0.5 + 
   sigma2^0.5* (-0.5 + 0.5* (exp(-eta2)* y2)^sqrt(sigma2))* log(
     exp(-eta2)* y2)))/sigma2^0.5
            

}






}

####



####

if(margin2 == "WEI2"){


          
  pdf2 <- exp(eta2)*sqrt(sigma2)*(exp(eta2)*y2)^(sqrt(sigma2)-1)*exp(-(exp(eta2)*y2)^sqrt(sigma2))

  
  
derpdf2.dereta2 <-  (exp(eta2) * sqrt(sigma2) * (exp(eta2) * y2)^(sqrt(sigma2) - 
    1) + exp(eta2) * sqrt(sigma2) * ((exp(eta2) * y2)^((sqrt(sigma2) - 
    1) - 1) * ((sqrt(sigma2) - 1) * (exp(eta2) * y2)))) * exp(-(exp(eta2) * 
    y2)^sqrt(sigma2)) - exp(eta2) * sqrt(sigma2) * (exp(eta2) * 
    y2)^(sqrt(sigma2) - 1) * (exp(-(exp(eta2) * y2)^sqrt(sigma2)) * 
    ((exp(eta2) * y2)^(sqrt(sigma2) - 1) * (sqrt(sigma2) * (exp(eta2) * 
        y2))))
                    
derpdf2.sigma2 <- (exp(eta2) * (0.5 * sigma2^-0.5) * (exp(eta2) * y2)^(sqrt(sigma2) - 
    1) + exp(eta2) * sqrt(sigma2) * ((exp(eta2) * y2)^(sqrt(sigma2) - 
    1) * (log((exp(eta2) * y2)) * (0.5 * sigma2^-0.5)))) * exp(-(exp(eta2) * 
    y2)^sqrt(sigma2)) - exp(eta2) * sqrt(sigma2) * (exp(eta2) * 
    y2)^(sqrt(sigma2) - 1) * (exp(-(exp(eta2) * y2)^sqrt(sigma2)) * 
    ((exp(eta2) * y2)^sqrt(sigma2) * (log((exp(eta2) * y2)) * 
        (0.5 * sigma2^-0.5))))


dersigma2.dersigma2.st <- exp(sigma2.st)   

 der2pdf2.dereta2 <- (sqrt(sigma2)*(sigma2*exp(3*eta2*sqrt(sigma2))*y2^(3*sqrt(sigma2))-3*sigma2*exp(2*eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))+sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*exp(-exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)))/y2 
    
    
    
der2pdf2.dersigma22 <- -(exp(-exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*((-sigma2*exp(3*eta2*sqrt(sigma2))*y2^(3*sqrt(sigma2))+3*sigma2*exp(2*eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))-sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*
log(exp(eta2)*y2)^2+sqrt(sigma2)*(exp(2*eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))-exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*log(exp(eta2)*y2)+exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2)))/(4*sigma2^(3/2)*y2)    
    
    
    
    
der2pdf2.dereta2dersigma2 <- -(exp(-exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*((-sigma2*exp(3*eta2*sqrt(sigma2))*y2^(3*sqrt(sigma2))+3*sigma2*exp(2*eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))-sigma2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))*
log(exp(eta2)*y2)+sqrt(sigma2)*(2*exp(2*eta2*sqrt(sigma2))*y2^(2*sqrt(sigma2))-2*exp(eta2*sqrt(sigma2))*y2^sqrt(sigma2))))/(2*sqrt(sigma2)*y2)
 
 
 
 
if(naive == FALSE){   
 
    
    p2  <-  1  - exp( - (exp(eta2)*y2)^sqrt(sigma2) )


derp2.dereta2 <- exp(-(exp(eta2) * y2)^sqrt(sigma2)) * ((exp(eta2) * y2)^(sqrt(sigma2) - 
    1) * (sqrt(sigma2) * (exp(eta2) * y2)))

                
                          
derp2.dersigma2 <- exp(-(exp(eta2) * y2)^sqrt(sigma2)) * ((exp(eta2) * y2)^sqrt(sigma2) * 
    (log((exp(eta2) * y2)) * (0.5 * sigma2^-0.5)))           
    


der2p2.dereta2eta2 <- exp(-(exp(eta2) * y2)^sqrt(sigma2)) * ((exp(eta2) * y2)^((sqrt(sigma2) - 
    1) - 1) * ((sqrt(sigma2) - 1) * (exp(eta2) * y2)) * (sqrt(sigma2) * 
    (exp(eta2) * y2)) + (exp(eta2) * y2)^(sqrt(sigma2) - 1) * 
    (sqrt(sigma2) * (exp(eta2) * y2))) - exp(-(exp(eta2) * y2)^sqrt(sigma2)) * 
    ((exp(eta2) * y2)^(sqrt(sigma2) - 1) * (sqrt(sigma2) * (exp(eta2) * 
        y2))) * ((exp(eta2) * y2)^(sqrt(sigma2) - 1) * (sqrt(sigma2) * 
    (exp(eta2) * y2)))



der2p2.dersigma22 <- exp(-(exp(eta2) * y2)^sqrt(sigma2)) * ((exp(eta2) * y2)^sqrt(sigma2) * 
    (log((exp(eta2) * y2)) * (0.5 * sigma2^-0.5)) * (log((exp(eta2) * 
    y2)) * (0.5 * sigma2^-0.5)) + (exp(eta2) * y2)^sqrt(sigma2) * 
    (log((exp(eta2) * y2)) * (0.5 * (-0.5 * sigma2^-1.5)))) - 
    exp(-(exp(eta2) * y2)^sqrt(sigma2)) * ((exp(eta2) * y2)^sqrt(sigma2) * 
        (log((exp(eta2) * y2)) * (0.5 * sigma2^-0.5))) * ((exp(eta2) * 
        y2)^sqrt(sigma2) * (log((exp(eta2) * y2)) * (0.5 * sigma2^-0.5)))


der2p2.dereta2dersigma2 <- exp(-(exp(eta2) * y2)^sqrt(sigma2)) * ((exp(eta2) * y2)^(sqrt(sigma2) - 
    1) * (sqrt(sigma2) * (exp(eta2) * y2)) * (log((exp(eta2) * 
    y2)) * (0.5 * sigma2^-0.5)) + (exp(eta2) * y2)^sqrt(sigma2) * 
    (exp(eta2) * y2/(exp(eta2) * y2) * (0.5 * sigma2^-0.5))) - 
    exp(-(exp(eta2) * y2)^sqrt(sigma2)) * ((exp(eta2) * y2)^(sqrt(sigma2) - 
        1) * (sqrt(sigma2) * (exp(eta2) * y2))) * ((exp(eta2) * 
        y2)^sqrt(sigma2) * (log((exp(eta2) * y2)) * (0.5 * sigma2^-0.5)))

       
}






}

####









if(margin2 == "iG"){

  pdf2          <- exp(-0.5 * log(2 * pi) - log(sqrt(sigma2)) - (3/2) * log(y2) - 
                   ((y2 - exp(eta2))^2)/(2 * sigma2 * (exp(eta2)^2) * y2))
                  
                  
 derpdf2.dereta2 <-  exp(-0.5 * log(2 * pi) - log(sqrt(sigma2)) - (3/2) * log(y2) - 
                    ((y2 - exp(eta2))^2)/(2 * sigma2 * (exp(eta2)^2) * y2)) * 
                    (2 * (exp(eta2) * (y2 - exp(eta2)))/(2 * sigma2 * (exp(eta2)^2) * 
                    y2) + ((y2 - exp(eta2))^2) * (2 * sigma2 * (2 * (exp(eta2) * 
                    exp(eta2))) * y2)/(2 * sigma2 * (exp(eta2)^2) * y2)^2)                 
                  
                  
   derpdf2.sigma2 <- -(exp(-0.5 * log(2 * pi) - log(sqrt(sigma2)) - (3/2) * log(y2) - 
                   ((y2 - exp(eta2))^2)/(2 * sigma2 * (exp(eta2)^2) * y2)) * 
                   (0.5 * sigma2^-0.5/sqrt(sigma2) - ((y2 - exp(eta2))^2) * 
                   (2 * (exp(eta2)^2) * y2)/(2 * sigma2 * (exp(eta2)^2) * 
                    y2)^2))


dersigma2.dersigma2.st <- exp(sigma2.st)                
                  
 der2pdf2.dereta2 <-  ((exp(exp(-eta2)/sigma2)*y2^2+(-2*exp(2*eta2)*sigma2-2*exp(eta2))*exp(exp(-eta2)/sigma2)*y2+(exp(3*eta2)*sigma2+exp(2*eta2))*exp(exp(-eta2)/sigma2))*
exp(-(3*log(y2))/2-(exp(-2*eta2)*y2)/(2*sigma2)-1/(2*sigma2*y2)-log(sigma2)/2-log(2*pi)/2-4*eta2))/sigma2^2 
   
   
der2pdf2.dersigma22 <- ((exp(exp(-eta2)/sigma2)*y2^4+(-6*exp(2*eta2)*sigma2-4*exp(eta2))*exp(exp(-eta2)/sigma2)*y2^3+(3*exp(4*eta2)*sigma2^2+12*exp(3*eta2)*sigma2+6*exp(2*eta2))*exp(exp(-eta2)/sigma2)*y2^2+
(-6*exp(4*eta2)*sigma2-4*exp(3*eta2))*exp(exp(-eta2)/sigma2)*y2+exp(exp(-eta2)/sigma2+4*eta2))*exp(-(3*log(y2))/2-(exp(-2*eta2)*y2)/(2*sigma2)-1/(2*sigma2*y2)-log(sigma2)/2-log(2*pi)/2-4*eta2))/(4*sigma2^4*y2^2)
 
 
der2pdf2.dereta2dersigma2 <-  ((exp(exp(-eta2)/sigma2)*y2^3+(-3*exp(2*eta2)*sigma2-3*exp(eta2))*exp(exp(-eta2)/sigma2)*y2^2+(3*exp(3*eta2)*sigma2+3*exp(2*eta2))*exp(exp(-eta2)/sigma2)*y2-exp(exp(-eta2)/sigma2+3*eta2))*
exp(-(3*log(y2))/2-(exp(-2*eta2)*y2)/(2*sigma2)-1/(2*sigma2*y2)-log(sigma2)/2-log(2*pi)/2-4*eta2))/(2*sigma2^3*y2)                
                  
                  
if(naive == FALSE){                   
          
                  
    p2          <-  pnorm(((y2/exp(eta2)) - 1)/(sqrt(sigma2) * sqrt(y2))) + 
                    exp(2/(exp(eta2)*sigma2))* pnorm(-((y2/exp(eta2)) + 1)/(sqrt(sigma2) * sqrt(y2)))
                



derp2.dereta2    <- exp(2/(exp(eta2) * sigma2)) * (dnorm(-((y2/exp(eta2)) + 1)/(sqrt(sigma2) * 
                    sqrt(y2))) * (y2 * exp(eta2)/exp(eta2)^2/(sqrt(sigma2) * 
                    sqrt(y2)))) - exp(2/(exp(eta2) * sigma2)) * (2 * (exp(eta2) * 
                    sigma2)/(exp(eta2) * sigma2)^2) * pnorm(-((y2/exp(eta2)) + 
                    1)/(sqrt(sigma2) * sqrt(y2))) - dnorm(((y2/exp(eta2)) - 1)/(sqrt(sigma2) * 
                    sqrt(y2))) * (y2 * exp(eta2)/exp(eta2)^2/(sqrt(sigma2) * 
                    sqrt(y2)))

                     


                          
derp2.dersigma2 <- exp(2/(exp(eta2) * sigma2)) * (dnorm(-((y2/exp(eta2)) + 1)/(sqrt(sigma2) * 
                   sqrt(y2))) * (((y2/exp(eta2)) + 1) * (0.5 * sigma2^-0.5 * 
                   sqrt(y2))/(sqrt(sigma2) * sqrt(y2))^2)) - exp(2/(exp(eta2) * 
                   sigma2)) * (2 * exp(eta2)/(exp(eta2) * sigma2)^2) * pnorm(-((y2/exp(eta2)) + 
                   1)/(sqrt(sigma2) * sqrt(y2))) - dnorm(((y2/exp(eta2)) - 1)/(sqrt(sigma2) * 
                   sqrt(y2))) * (((y2/exp(eta2)) - 1) * (0.5 * sigma2^-0.5 * 
                   sqrt(y2))/(sqrt(sigma2) * sqrt(y2))^2)                    
    


der2p2.dereta2eta2 <- (exp(-3*eta2)*(sqrt(sigma2)*(2*exp(2*eta2)*sigma2+4*exp(eta2))*exp((2*exp(-eta2))/sigma2)*pnorm(-(exp(-eta2)*(y2+exp(eta2)))/(sqrt(sigma2)*sqrt(y2)))+sqrt(y2)*(
(sigma2*exp((2*exp(-eta2))/sigma2)*y2+(-exp(2*eta2)*sigma2^2-3*exp(eta2)*sigma2)*exp((2*exp(-eta2))/sigma2))*dnorm(-(exp(-eta2)*(y2+exp(eta2)))/(sqrt(sigma2)*sqrt(y2)))+
(-sigma2*y2+exp(2*eta2)*sigma2^2+exp(eta2)*sigma2)*dnorm((exp(-eta2)*(y2-exp(eta2)))/(sqrt(sigma2)*sqrt(y2))))))/sigma2^(5/2)

    
der2p2.dersigma22 <-  (exp(-3*eta2)*(sqrt(sigma2)*(16*exp(2*eta2)*sigma2+16*exp(eta2))*exp((2*exp(-eta2))/sigma2)*y2^(3/2)*pnorm(-(exp(-eta2)*(y2+exp(eta2)))/(sqrt(sigma2)*sqrt(y2)))+(sigma2*exp((2*exp(-eta2))/sigma2)*y2^3+
(-3*exp(2*eta2)*sigma2^2-5*exp(eta2)*sigma2)*exp((2*exp(-eta2))/sigma2)*y2^2+(-3*exp(3*eta2)*sigma2^2-5*exp(2*eta2)*sigma2)*exp((2*exp(-eta2))/sigma2)*y2+sigma2*exp((2*exp(-eta2))/sigma2+3*eta2))*
dnorm(-(exp(-eta2)*(y2+exp(eta2)))/(sqrt(sigma2)*sqrt(y2)))+
(-sigma2*y2^3+(3*exp(2*eta2)*sigma2^2+3*exp(eta2)*sigma2)*y2^2+(-3*exp(3*eta2)*sigma2^2-3*exp(2*eta2)*sigma2)*y2+exp(3*eta2)*sigma2)*dnorm((exp(-eta2)*(y2-exp(eta2)))/(sqrt(sigma2)*sqrt(y2)))
))/(4*sigma2^(9/2)*y2^(3/2))




der2p2.dereta2dersigma2 <- (exp(-3*eta2)*((4*exp(2*eta2)*sigma2+8*exp(eta2))*exp((2*exp(-eta2))/sigma2)*sqrt(y2)*pnorm(-(exp(-eta2)*(y2+exp(eta2)))/(sqrt(sigma2)*sqrt(y2)))+sqrt(sigma2)*(
(exp((2*exp(-eta2))/sigma2)*y2^2+(-exp(2*eta2)*sigma2-4*exp(eta2))*exp((2*exp(-eta2))/sigma2)*y2-exp((2*exp(-eta2))/sigma2+2*eta2))*dnorm(-(exp(-eta2)*(y2+exp(eta2)))/(sqrt(sigma2)*sqrt(y2)))+
(-y2^2+(exp(2*eta2)*sigma2+2*exp(eta2))*y2-exp(2*eta2))*dnorm((exp(-eta2)*(y2-exp(eta2)))/(sqrt(sigma2)*sqrt(y2))))))/(2*sigma2^3*sqrt(y2))



}





}


####

if(margin2 == "LO"){

  pdf2          <- dlogis(y2,eta2,sqrt(sigma2)) # exp(-(y2-eta2)/sqrt(sigma2))/(sqrt(sigma2)*(1+exp(-(y2-eta2)/sqrt(sigma2)))^2)
  
derpdf2.dereta2 <-     (exp((eta2 + y2)/sqrt(sigma2))* (-exp((eta2/sqrt(sigma2))) + exp(y2/sqrt(
                        sigma2))))/((exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(sigma2)))^3 *sigma2)  
  
derpdf2.sigma2 <-             
     (exp((eta2 + y2)/sqrt(
   sigma2))* (exp(eta2/sqrt(
      sigma2))* (0.5* eta2* sigma2^2.5 - 0.5* sigma2^3. - 
        0.5 *sigma2^2.5* y2) + 
     exp(y2/sqrt(
      sigma2))* (-0.5* eta2 *sigma2^2.5 - 0.5 *sigma2^3 + 
        0.5* sigma2^2.5* y2)))/((exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(
     sigma2)))^3 *sigma2^4.5)  
  
dersigma2.dersigma2.st <- exp(sigma2.st)   
  
 der2pdf2.dereta2 <- (exp((eta2 + y2)/sqrt(sigma2))* (exp((2* eta2)/sqrt(sigma2)) + exp((2 *y2)/sqrt(
   sigma2)) - 4* exp((eta2 + y2)/sqrt(sigma2))))/((exp(eta2/sqrt(
   sigma2)) + exp(y2/sqrt(sigma2)))^4* sigma2^(3/2))  
  
der2pdf2.dersigma22 <- (1/((exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(sigma2)))^4 *sigma2^13.5))*exp((
 eta2 + y2)/sqrt(sigma2))* (exp((eta2 + y2)/sqrt(
    sigma2))* (-1* eta2^2* sigma2^10 + 1.5 *sigma2^11 + 
      2* eta2* sigma2^10* y2 - 1*sigma2^10* y2^2) + 
   exp((2* y2)/sqrt(
    sigma2))* (0.25 *eta2^2* sigma2^10 + 1.25* eta2* sigma2^10.5 + 
      0.75 *sigma2^11 - 0.5* eta2* sigma2^10* y2 - 
      1.25* sigma2^10.5* y2 + 0.25* sigma2^10* y2^2) + 
   exp((2* eta2)/sqrt(
    sigma2))* (0.25* eta2^2* sigma2^10 - 1.25* eta2* sigma2^10.5 + 
      0.75* sigma2^11 - 0.5* eta2* sigma2^10* y2 + 
      1.25* sigma2^10.5* y2 + 0.25* sigma2^10* y2^2))  
  
der2pdf2.dereta2dersigma2 <- (exp((eta2 + 
  y2)/sqrt(sigma2))* (exp((eta2 + y2)/sqrt(sigma2))*
     sigma2^4* (2* eta2 - 2* y2) + 
   exp((2* y2)/sqrt(
    sigma2))* (-0.5* eta2 *sigma2^4 - 1* sigma2^4.5 + 
      0.5* sigma2^4* y2) + 
   exp((2 *eta2)/sqrt(
    sigma2))* (-0.5* eta2* sigma2^4 + 1* sigma2^4.5 + 
      0.5* sigma2^4* y2)))/((exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(
   sigma2)))^4* sigma2^6.5)  
  
  
  
  
  
  
  
  
  
  
  
  
if(naive == FALSE){   
  
    p2          <- plogis(y2,eta2,sqrt(sigma2)) #1/(1+exp(-(y2-eta2)/sqrt(sigma2)))
             



derp2.dereta2    <- -(exp(-(y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))/(1 + exp(-(y2 - 
                      eta2)/sqrt(sigma2)))^2)
                    




                          
derp2.dersigma2 <-  -(exp(-(y2 - eta2)/sqrt(sigma2)) * ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)/(1 + 
                      exp(-(y2 - eta2)/sqrt(sigma2)))^2)                     
    


der2p2.dereta2eta2 <-  (exp((2* eta2 + y2)/sqrt(sigma2)) - exp((
 eta2 + 2* y2)/sqrt(sigma2)))/((exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(
   sigma2)))^3* sigma2)




   
der2p2.dersigma22 <-  (0.5 *exp((2* (eta2 - y2))/sqrt(
  sigma2))* (eta2 - y2)^2)/((1 + exp((eta2 - y2)/sqrt(
    sigma2)))^3* sigma2^3) - (
 0.25 * exp((eta2 + y2)/sqrt(
  sigma2)) *(eta2 - 1*y2)* (1*eta2* sigma2^2.5 + 3* sigma2^3 - 
    1* sigma2^2.5* y2))/((exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(
    sigma2)))^2* sigma2^5.5)






der2p2.dereta2dersigma2 <- (exp((eta2 + 
  y2)/sqrt(sigma2))* (exp(y2/sqrt(
    sigma2))* (0.5* eta2* sigma2^1.5 + 0.5* sigma2^2 - 
      0.5* sigma2^1.5 *y2) + 
   exp(eta2/sqrt(
    sigma2))* (-0.5 *eta2* sigma2^1.5 + 0.5 *sigma2^2 + 
      0.5* sigma2^1.5 *y2)))/((exp(eta2/sqrt(sigma2)) + exp(y2/sqrt(
   sigma2)))^3 *sigma2^3.5)
            

}



}









if(margin2 == "rGU"){

#sigma2    <- ifelse(sigma2 < 0.006, 0.006, sigma2)
#sigma2.st <- log(sigma2) 

  pdf2          <- 1/sqrt(sigma2)*exp(-((y2-eta2)/sqrt(sigma2)+exp(-((y2-eta2)/sqrt(sigma2)))))
  
derpdf2.dereta2 <-  -(1/sqrt(sigma2) * (exp(-((y2 - eta2)/sqrt(sigma2) + exp(-((y2 - 
                      eta2)/sqrt(sigma2))))) * (exp(-((y2 - eta2)/sqrt(sigma2))) * 
                      (1/sqrt(sigma2)) - 1/sqrt(sigma2))))  
  
derpdf2.sigma2 <-  (exp(-exp(((eta2 - y2)/sqrt(sigma2))) + (eta2 - 2 * y2)/sqrt(
  sigma2)) * (exp(eta2/sqrt(sigma2))* (0.5 * eta2 - 0.5 * y2) + 
   exp(y2/sqrt(sigma2))*(-0.5* eta2 - 0.5* sqrt(sigma2) + 0.5 * y2)))/sigma2^2

dersigma2.dersigma2.st <- exp(sigma2.st)   
  
 der2pdf2.dereta2 <-  -(1/sqrt(sigma2) * (exp(-((y2 - eta2)/sqrt(sigma2) + exp(-((y2 - 
                       eta2)/sqrt(sigma2))))) * (exp(-((y2 - eta2)/sqrt(sigma2))) * 
                       (1/sqrt(sigma2)) * (1/sqrt(sigma2))) - exp(-((y2 - eta2)/sqrt(sigma2) + 
                        exp(-((y2 - eta2)/sqrt(sigma2))))) * (exp(-((y2 - eta2)/sqrt(sigma2))) * 
                       (1/sqrt(sigma2)) - 1/sqrt(sigma2)) * (exp(-((y2 - eta2)/sqrt(sigma2))) * 
                       (1/sqrt(sigma2)) - 1/sqrt(sigma2))))  
  
der2pdf2.dersigma22 <-  (1/(sigma2^9))*exp(-exp(((eta2 - y2)/sqrt(sigma2))) + (eta2 - 3* y2)/
  sqrt(sigma2))* (exp((2 * eta2)/sqrt(sigma2))*
     sigma2^5.5* (0.25* eta2^2 - 0.5* eta2* y2 + 0.25* y2^2) + 
   exp((eta2 + y2)/sqrt(
    sigma2))* (-0.75* eta2^2 *sigma2^5.5 - 1.25* eta2* sigma2^6 + 
      1.5* eta2* sigma2^5.5* y2 + 1.25* sigma2^6* y2 - 
      0.75* sigma2^5.5* y2^2) + 
   exp((2 *y2)/sqrt(
    sigma2))* (0.25* eta2^2 *sigma2^5.5 + 1.25* eta2* sigma2^6 + 
      0.75* sigma2^6.5 - 0.5* eta2* sigma2^5.5* y2 - 1.25* sigma2^6* y2 + 
      0.25* sigma2^5.5* y2^2))  
  
der2pdf2.dereta2dersigma2 <- (1/(sigma2^4.5))*exp(-exp(((eta2 - y2)/sqrt(sigma2))) + (eta2 - 3 *y2)/
  sqrt(sigma2)) *(exp((2* eta2)/sqrt(sigma2))*
     sigma2^2* (-0.5* eta2 + 0.5 *y2) + 
   exp((eta2 + y2)/sqrt(
    sigma2))* (1.5* eta2* sigma2^2 + sigma2^2.5 - 1.5* sigma2^2* y2) +
    exp((2* y2)/sqrt(
    sigma2))* (-0.5* eta2* sigma2^2 - sigma2^2.5 + 0.5* sigma2^2* y2))  
  
  
  
  
  
  
  
if(naive == FALSE){   
  
  
  
    p2          <- exp(-(exp(-(y2-eta2)/sqrt(sigma2))))
                


derp2.dereta2    <- -(exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                      (1/sqrt(sigma2))))
                      
                      
                   


                          
derp2.dersigma2 <-   -(exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                     ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)))                      
    


der2p2.dereta2eta2 <- -(exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                       (1/sqrt(sigma2)) * (1/sqrt(sigma2))) - exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * 
                       (exp(-(y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))) * (exp(-(y2 - 
                        eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))))




der2p2.dersigma22 <- exp(-exp(((eta2 - y2)/sqrt(sigma2))) + (eta2 - 2 *y2)/sqrt(
  sigma2)) *(eta2 - 1* y2)* ((
   0.25* exp(eta2/sqrt(sigma2))* (eta2 - y2))/sigma2^3 + 
   exp(y2/sqrt(sigma2))* (-((0.25 * eta2)/sigma2^3) - 0.75/sigma2^2.5 + (0.25 * y2)/
      sigma2^3))






der2p2.dereta2dersigma2 <- -(exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                            ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2) * (1/sqrt(sigma2)) - 
                             exp(-(y2 - eta2)/sqrt(sigma2)) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2)) - 
                             exp(-(exp(-(y2 - eta2)/sqrt(sigma2)))) * (exp(-(y2 - eta2)/sqrt(sigma2)) * 
                            ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)) * 
                            (exp(-(y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))))
            

}




}




if(margin2 == "GU"){



  pdf2          <- exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2)))
                   
 derpdf2.dereta2 <-  exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                    (1/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))) - 
                     exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                      (1/sqrt(sigma2)) * (1/sqrt(sigma2)))                   
                   
derpdf2.sigma2 <-  (exp(-exp(((-eta2 + y2)/sqrt(sigma2))) + (-2* eta2 + y2)/sqrt(
  sigma2)) *(exp(y2/sqrt(sigma2)) *sigma2^1.5* (-0.5* eta2 + 0.5* y2) + 
   exp(eta2/sqrt(
    sigma2))* (0.5* eta2* sigma2^1.5 - 0.5 *sigma2^2 - 
      0.5* sigma2^1.5* y2)))/sigma2^3.5 

dersigma2.dersigma2.st <- exp(sigma2.st)                    
                   
der2pdf2.dereta2 <-  (exp(-exp(((-eta2 + y2)/sqrt(sigma2))) + (-3 *eta2 + y2)/sqrt(
  sigma2)) *(exp((2 *eta2)/sqrt(sigma2)) + exp((2* y2)/sqrt(sigma2)) - 
   3* exp((eta2 + y2)/sqrt(sigma2))))/sigma2^(3/2)                   
                   
der2pdf2.dersigma22 <- (1/(sigma2^11))*exp(-exp(((-eta2 + y2)/sqrt(sigma2))) + (-3 *eta2 + y2)/
  sqrt(sigma2))* (exp((2* y2)/sqrt(sigma2))*
     sigma2^7.5* (0.25* eta2^2 - 0.5 *eta2* y2 + 0.25* y2^2) + 
   exp((eta2 + y2)/sqrt(
    sigma2))* (-0.75 *eta2^2* sigma2^7.5 + 1.25* eta2* sigma2^8 + 
      1.5 *eta2* sigma2^7.5* y2 - 1.25* sigma2^8* y2 - 
      0.75* sigma2^7.5* y2^2) + 
   exp((2* eta2)/sqrt(
    sigma2))* (0.25* eta2^2* sigma2^7.5 - 1.25 *eta2* sigma2^8 + 
      0.75* sigma2^8.5 - 0.5* eta2* sigma2^7.5* y2 + 1.25* sigma2^8* y2 + 
      0.25* sigma2^7.5* y2^2))                   
                   
der2pdf2.dereta2dersigma2 <-(1/(sigma2^6.5))*exp(-exp(((-eta2 + y2)/sqrt(sigma2))) + (-3* eta2 + y2)/
  sqrt(sigma2)) *(exp((2* y2)/sqrt(sigma2))*
     sigma2^4* (-0.5* eta2 + 0.5* y2) + 
   exp((eta2 + y2)/sqrt(
    sigma2)) *(1.5* eta2* sigma2^4 - sigma2^4.5 - 1.5* sigma2^4* y2) +
    exp((2* eta2)/sqrt(
    sigma2))* (-0.5* eta2* sigma2^4 +  sigma2^4.5 + 0.5* sigma2^4* y2))                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
if(naive == FALSE){                    
                   
    p2          <- 1 - exp(-exp((y2 - eta2)/sqrt(sigma2)))
    
 

derp2.dereta2    <- -(exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                     (1/sqrt(sigma2))))
                  


                          
derp2.dersigma2 <-   -(exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                      ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)))                       
    


der2p2.dereta2eta2 <- -(exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                        (1/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))) - 
                        exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                        (1/sqrt(sigma2)) * (1/sqrt(sigma2))))



 

 der2p2.dersigma22 <-  exp(-exp(((-eta2 + y2)/sqrt(sigma2))) + (-2* eta2 + y2)/sqrt(
  sigma2))* (-((0.25* exp(y2/sqrt(sigma2))* (-eta2 + y2)^2)/
    sigma2^3) + (
   0.25 *exp(eta2/sqrt(
    sigma2)) *(eta2 - y2)* (eta2* sigma2^2.5 - 3* sigma2^3 - 
      sigma2^2.5* y2))/sigma2^5.5)






der2p2.dereta2dersigma2 <- -(exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                            ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)) * (exp((y2 - 
                              eta2)/sqrt(sigma2)) * (1/sqrt(sigma2))) - exp(-exp((y2 - 
                              eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * (0.5 * 
                              sigma2^-0.5/sqrt(sigma2)^2) + exp((y2 - eta2)/sqrt(sigma2)) * 
                            ((y2 - eta2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2) * (1/sqrt(sigma2))))
            

}




}



if(margin2 %in% c("GA")){


sigma2    <- ifelse(sigma2 < 0.006, 0.006, sigma2)
sigma2.st <- log(sigma2) 


 pdf2  <-  dgamma(y2, shape = 1/sigma2, scale = exp(eta2) * sigma2)
           
  derpdf2.dereta2 <- -((exp(eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1 - 1/sigma2)*
    y2^(-1 + 1/sigma2))/gamma(1/sigma2)) + (
 exp(-eta2 - (exp(-eta2)* y2)/sigma2) * (exp(eta2)* sigma2)^(-1/sigma2)* y2^(1/
  sigma2))/(sigma2* gamma(1/sigma2))         
           
derpdf2.sigma2 <- (exp(-eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1/sigma2)* y2^(1/
  sigma2))/(sigma2^2* gamma(1/sigma2)) + (
 exp(-((exp(-eta2)* y2)/sigma2))* (exp(eta2)* sigma2)^(-1/sigma2)*
   y2^(-1 + 1/sigma2)* (-(1/sigma2^2) + log(exp(eta2)* sigma2)/sigma2^2))/
 gamma(1/sigma2) - (
 exp(-((exp(-eta2)* y2)/sigma2))* (exp(eta2)* sigma2)^(-1/sigma2)*
   y2^(-1 + 1/sigma2)* log(y2))/(sigma2^2* gamma(1/sigma2)) + (
 exp(-((exp(-eta2)* y2)/sigma2))* (exp(eta2)* sigma2)^(-1/sigma2)*
   y2^(-1 + 1/sigma2)* psigamma(1/sigma2,0))/(
 sigma2^2 *gamma(1/sigma2))
 
 

dersigma2.dersigma2.st <- exp(sigma2.st)           
           

 der2pdf2.dereta2 <-  (y2^(-1 + 1/
  sigma2)* (exp(-((exp(-eta2)* y2)/
     sigma2))* (-exp(
        2* eta2) *(-1 - 1/sigma2)* sigma2 *(exp(eta2)* sigma2)^(-2 - 1/
        sigma2) - exp(eta2)* (exp(eta2)* sigma2)^(-1 - 1/sigma2)) - (
   2* exp(-((exp(-eta2)* y2)/sigma2))* (exp(eta2)* sigma2)^(-1 - 1/sigma2)* y2)/
   sigma2 + (
   exp(-eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1/sigma2)*
     y2* (-1 + (exp(-eta2)* y2)/sigma2))/sigma2))/gamma(1/sigma2)           
           
           
           
der2pdf2.dersigma22 <- (1/gamma(1/sigma2))*
 y2^(-1 + 1/
   sigma2)* ((exp(eta2)* sigma2)^(-1/
      sigma2)* (-((2* exp(-eta2 - (exp(-eta2)* y2)/sigma2)* y2)/sigma2^3) + (
       exp(-2* eta2 - (exp(-eta2)* y2)/sigma2)* y2^2)/sigma2^4) + (
    2 * exp(-eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1/sigma2)*
      y2* (-(1/sigma2^2) + log(exp(eta2)* sigma2)/sigma2^2))/sigma2^2 + 
    exp(-((exp(-eta2)* y2)/
      sigma2))* ((exp(eta2)* sigma2)^(-1/
         sigma2)* (3/sigma2^3 - (2* log(exp(eta2)* sigma2))/sigma2^3) + (exp(
          eta2)* sigma2)^(-1/
         sigma2)* (-(1/sigma2^2) + log(exp(eta2)* sigma2)/sigma2^2)^2)) + 
 2 *((exp(-eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1/sigma2)* y2)/
    sigma2^2 + 
    exp(-((exp(-eta2)* y2)/sigma2))* (exp(eta2)* sigma2)^(-1/
      sigma2)* (-(1/sigma2^2) + log(exp(eta2)* sigma2)/sigma2^2))* (-((
     y2^(-1 + 1/sigma2)* log(y2))/(sigma2^2 *gamma(1/sigma2))) + (
    y2^(-1 + 1/sigma2)* psigamma(1/sigma2,0))/(
    sigma2^2 *gamma(1/sigma2))) + 
 exp(-((exp(-eta2)* y2)/sigma2))* (exp(eta2)* sigma2)^(-1/
   sigma2)* (((2* y2^(-1 + 1/sigma2)* log(y2))/sigma2^3 + (
     y2^(-1 + 1/sigma2)* log(y2)^2)/sigma2^4)/gamma(1/sigma2) - (
    2 *y2^(-1 + 1/sigma2)* log(y2)* psigamma(1/sigma2,0))/(
    sigma2^4* gamma(1/sigma2)) + 
    y2^(-1 + 1/
      sigma2)* (-((2 *psigamma(1/sigma2,0))/(
        sigma2^3* gamma(1/sigma2))) + psigamma(1/sigma2,0)^2/(
       sigma2^4 *gamma(1/sigma2)) - psigamma(1/sigma2,1)/(
       sigma2^4* gamma(1/sigma2))))           
           
           
der2pdf2.dereta2dersigma2 <- (exp(-eta2 - (exp(-eta2)* y2)/sigma2) * (exp(eta2)* sigma2)^(-1/
   sigma2) * (exp(eta2) - y2) * y2^(-1 + 1/sigma2))/(
 sigma2^2 * gamma(1/sigma2)) - (
 exp(-2* eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1/
   sigma2)* (exp(eta2) - y2)* y2^(1/sigma2))/(sigma2^3 * gamma(1/sigma2)) - (
 exp(-eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1/
   sigma2)* (exp(eta2) - y2)* y2^(-1 + 1/
   sigma2)* (-(1/sigma2^2) + log(exp(eta2)* sigma2)/sigma2^2))/(
 sigma2* gamma(1/sigma2)) + (
 exp(-eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1/
   sigma2)* (exp(eta2) - y2)* y2^(-1 + 1/sigma2)* log(y2))/(
 sigma2^3* gamma(1/sigma2)) - (
 exp(-eta2 - (exp(-eta2)* y2)/sigma2)* (exp(eta2)* sigma2)^(-1/
   sigma2)* (exp(eta2) - y2)* y2^(-1 + 1/sigma2)* psigamma(1/sigma2,0))/(
 sigma2^3 *gamma(1/sigma2))           
           
           
    
    
    
if(naive == FALSE){      
    
           
    p2          <-  pgamma(y2, shape = 1/sigma2, scale = exp(eta2) * sigma2)
                                  

derp2.dereta2    <- -((exp(-eta2 - (exp(-eta2)* y2)/sigma2)*
   y2* ((exp(-eta2)* y2)/sigma2)^(-1 + 1/sigma2))/(sigma2 *gamma(1/sigma2)))
   
   
der2p2.dereta2eta2 <- (exp(-2* eta2 - (exp(-eta2)* y2)/sigma2)* (-1 + 1/sigma2)* y2^2* ((exp(-eta2)* y2)/
   sigma2)^(-2 + 1/sigma2))/(sigma2^2* gamma(1/sigma2)) - (
 exp(-eta2 - (exp(-eta2)* y2)/sigma2)*
   y2* ((exp(-eta2)* y2)/sigma2)^(-1 + 1/
   sigma2)* (-1 + (exp(-eta2)* y2)/sigma2))/(sigma2* gamma(1/sigma2))
   

der2p2.dereta2dersigma2 <- -((exp(-eta2 - (exp(-eta2)* y2)/sigma2)* y2* ((exp(-eta2)* y2)/sigma2)^(1/
   sigma2))/(sigma2^2* gamma(1/sigma2))) - (
 exp(-((exp(-eta2)* y2)/sigma2))* ((exp(-eta2)* y2)/sigma2)^(1/
  sigma2)* (-(1/sigma2^2) - log((exp(-eta2)* y2)/sigma2)/sigma2^2))/
 gamma(1/sigma2) - (
 exp(-((exp(-eta2)* y2)/sigma2))* ((exp(-eta2)* y2)/sigma2)^(1/sigma2)*
   psigamma(1/sigma2,0))/(sigma2^2 *gamma(1/sigma2))
   
   
}



}





























if(margin2 == "iGA"){



pdf2          <-  exp(1/sigma2 * eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
                     (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 1))/y2))

derpdf2.dereta2 <- exp(1/sigma2 * eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
    (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 1))/y2)) * 
    (1/sigma2 - exp(eta2) * (1/sigma2 + 1)/y2)    
    
 derpdf2.sigma2 <- -(exp(1/sigma2 * eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
    (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 1))/y2)) * 
    (1/sigma2 * (1/sigma2^2/(1/sigma2 + 1)) + 1/sigma2^2 * log(1/sigma2 + 
        1) + 1/sigma2^2 * eta2 - 1/sigma2^2 * digamma(1/sigma2) - 
        1/sigma2^2 * log(y2) - exp(eta2) * (1/sigma2^2)/y2))   
    
 dersigma2.dersigma2.st <- exp(sigma2.st)    
    
  der2pdf2.dereta2 <-exp(1/sigma2 * eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
    (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 1))/y2)) * 
    (1/sigma2 - exp(eta2) * (1/sigma2 + 1)/y2) * (1/sigma2 - 
    exp(eta2) * (1/sigma2 + 1)/y2) - exp(1/sigma2 * eta2 + 1/sigma2 * 
    log(1/sigma2 + 1) - lgamma(1/sigma2) - (1/sigma2 + 1) * log(y2) - 
    ((exp(eta2) * (1/sigma2 + 1))/y2)) * (exp(eta2) * (1/sigma2 + 
    1)/y2)




der2pdf2.dersigma22 <- exp(1/sigma2 * eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
    (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 1))/y2)) * 
    (2 * sigma2/(sigma2^2)^2 * eta2 + (1/sigma2^2 * (1/sigma2^2/(1/sigma2 + 
        1)) + 2 * sigma2/(sigma2^2)^2 * log(1/sigma2 + 1) + (1/sigma2 * 
        (2 * sigma2/(sigma2^2)^2/(1/sigma2 + 1) - 1/sigma2^2 * 
            (1/sigma2^2)/(1/sigma2 + 1)^2) + 1/sigma2^2 * (1/sigma2^2/(1/sigma2 + 
        1)))) - (1/sigma2^2 * (1/sigma2^2 * trigamma(1/sigma2)) + 
        2 * sigma2/(sigma2^2)^2 * digamma(1/sigma2)) - 2 * sigma2/(sigma2^2)^2 * 
        log(y2) - exp(eta2) * (2 * sigma2/(sigma2^2)^2)/y2) + 
    exp(1/sigma2 * eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
        (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 
        1))/y2)) * (1/sigma2 * (1/sigma2^2/(1/sigma2 + 1)) + 
        1/sigma2^2 * log(1/sigma2 + 1) + 1/sigma2^2 * eta2 - 
        1/sigma2^2 * digamma(1/sigma2) - 1/sigma2^2 * log(y2) - 
        exp(eta2) * (1/sigma2^2)/y2) * (1/sigma2 * (1/sigma2^2/(1/sigma2 + 
        1)) + 1/sigma2^2 * log(1/sigma2 + 1) + 1/sigma2^2 * eta2 - 
        1/sigma2^2 * digamma(1/sigma2) - 1/sigma2^2 * log(y2) - 
        exp(eta2) * (1/sigma2^2)/y2)   
    
der2pdf2.dereta2dersigma2 <- -(exp(1/sigma2 * eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
    (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 1))/y2)) * 
    (1/sigma2^2 - exp(eta2) * (1/sigma2^2)/y2) + exp(1/sigma2 * 
    eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
    (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 1))/y2)) * 
    (1/sigma2 * (1/sigma2^2/(1/sigma2 + 1)) + 1/sigma2^2 * log(1/sigma2 + 
        1) + 1/sigma2^2 * eta2 - 1/sigma2^2 * digamma(1/sigma2) - 
        1/sigma2^2 * log(y2) - exp(eta2) * (1/sigma2^2)/y2) * 
    (1/sigma2 - exp(eta2) * (1/sigma2 + 1)/y2))     
    
    
    
    
 if(naive == FALSE){     
    
    
    p2          <-  1-pgamma(((exp(eta2) * (1/sigma2 + 1))/y2), shape = 1/sigma2, scale=1)
                    


derp2.dereta2    <- -(exp(eta2 - (exp(eta2)* (1 + 1/sigma2))/y2)* (1 + 1/sigma2)* ((
  exp(eta2)* (1 + 1/sigma2))/y2)^(-1 + 1/sigma2))/(y2 * gamma(1/sigma2))
   
   

der2p2.dereta2eta2 <- -((exp(2 *eta2 - (exp(eta2)* (1 + 1/sigma2))/
   y2)* (-1 + 1/sigma2) *(1 + 1/sigma2)^2 *((exp(eta2)* (1 + 1/sigma2))/
   y2)^(-2 + 1/sigma2))/(y2^2 *gamma(1/sigma2)) + (
 exp(eta2 - (exp(eta2)* (1 + 1/sigma2))/
   y2)* (1 + 1/sigma2)* (1 - (exp(eta2)* (1 + 1/sigma2))/y2)* ((
   exp(eta2)* (1 + 1/sigma2))/y2)^(-1 + 1/sigma2))/(y2* gamma(1/sigma2)))



der2p2.dereta2dersigma2 <- -((exp(eta2 - (exp(eta2)* (1 + 1/sigma2))/y2) *((exp(eta2)* (1 + 1/sigma2))/y2)^(
  1/sigma2))/(sigma2^2 *y2* gamma(1/sigma2)) + (
 exp(-((exp(eta2)* (1 + 1/sigma2))/y2))* ((exp(eta2)* (1 + 1/sigma2))/y2)^(1/
  sigma2)* (-(1/((1 + 1/sigma2)* sigma2^3)) - 
    log((exp(eta2)* (1 + 1/sigma2))/y2)/sigma2^2))/gamma(1/sigma2) + (
 exp(-((exp(eta2)* (1 + 1/sigma2))/y2))* ((exp(eta2)* (1 + 1/sigma2))/y2)^(1/
  sigma2)* psigamma(1/sigma2))/(sigma2^2* gamma(1/sigma2)))
   
   
            
  }






}





if(margin2 %in% c(cont2par,cont3par)){





derpdf2.dersigma2.st         <- derpdf2.sigma2 * dersigma2.dersigma2.st   
der2pdf2.dersigma2.st2       <- der2pdf2.dersigma22 * dersigma2.dersigma2.st^2 + derpdf2.sigma2  * dersigma2.dersigma2.st     
der2pdf2.dereta2dersigma2.st <- der2pdf2.dereta2dersigma2 *  dersigma2.dersigma2.st


if(margin2 %in% cont3par){

der2pdf2.dereta2dernu.st     <- der2pdf2.dereta2dernu * dernu.dernu.st
der2pdf2.sigma2.st2dernu.st  <- der2pdf2.dersigma2dernu * dersigma2.dersigma2.st * dernu.dernu.st                              
derpdf2.dernu.st             <- derpdf2.nu * dernu.dernu.st 
der2pdf2.dernu.st2           <- der2pdf2.dernu2 * dernu.dernu.st^2 +  derpdf2.nu  * dernu.dernu.st 

}



if(naive == FALSE){  

derp2.dersigma.st            <- derp2.dersigma2 *  dersigma2.dersigma2.st 
der2p2.dersigma2.st2         <- der2p2.dersigma22 * dersigma2.dersigma2.st^2 + derp2.dersigma2 * dersigma2.dersigma2.st
der2p2.dereta2dersigma2.st   <- der2p2.dereta2dersigma2 *  dersigma2.dersigma2.st  


if(margin2 %in% cont3par){

derp2.nu.st                  <- derp2.dernu *  dernu.dernu.st 
der2p2.dernu.st2             <- der2p2.dernu2 * dernu.dernu.st^2 + derp2.dernu * dernu.dernu.st
der2p2.dereta2dernu.st       <- der2p2.dereta2dernu * dernu.dernu.st 
der2p2.dersigma2.stdernu.st  <- der2p2.dersigma2dernu * dersigma2.dersigma2.st * dernu.dernu.st

}





}







}




epsilon <- 0.0000001 
max.p   <- 0.9999999

  pdf2 <- ifelse(pdf2 < sqrt(.Machine$double.eps), sqrt(.Machine$double.eps), pdf2 )

  p2   <- ifelse(p2 > max.p, max.p, p2) 
  p2   <- ifelse(p2 < epsilon, epsilon, p2) 



list(pdf2 = pdf2, p2 = p2, derpdf2.dereta2 = derpdf2.dereta2, 
     derpdf2.dersigma2.st = derpdf2.dersigma2.st, derp2.dersigma.st = derp2.dersigma.st,
     derp2.dereta2 = derp2.dereta2,
     der2p2.dereta2eta2 = der2p2.dereta2eta2, 
     der2pdf2.dereta2 = der2pdf2.dereta2,
     der2p2.dersigma2.st2 = der2p2.dersigma2.st2, 
     der2pdf2.dersigma2.st2 = der2pdf2.dersigma2.st2,
     der2p2.dereta2dersigma2.st = der2p2.dereta2dersigma2.st,            
     der2pdf2.dereta2dersigma2.st = der2pdf2.dereta2dersigma2.st,
     der2pdf2.dereta2dernu.st    = der2pdf2.dereta2dernu.st,   
     der2pdf2.sigma2.st2dernu.st = der2pdf2.sigma2.st2dernu.st,
     derpdf2.dernu.st            = derpdf2.dernu.st,           
     der2pdf2.dernu.st2          = der2pdf2.dernu.st2,         
     derp2.nu.st                 = derp2.nu.st,                
     der2p2.dernu.st2            = der2p2.dernu.st2,           
     der2p2.dereta2dernu.st      = der2p2.dereta2dernu.st,     
     der2p2.dersigma2.stdernu.st = der2p2.dersigma2.stdernu.st )     


}




    
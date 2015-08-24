copgHs <- function(p1, p2, eta1 = NULL, eta2 = NULL, teta, teta.st, BivD){



########################################################################################

if(!(BivD %in% c("N", "F")) ) {

derteta.derteta.st <- der2teta.derteta.stteta.st <- exp(teta.st) 
   
}   
   
########################################################################################
# Rotations
########################################################################################

if(BivD %in% c("C90","J90","G90") ) {
p1 <- 1 - p1 
teta <- -teta
}  

if(BivD %in% c("C180","J180","G180") ) {
p1 <- 1 - p1
p2 <- 1 - p2
}  

if(BivD %in% c("C270","J270","G270") ) {
p2 <- 1 - p2 
teta <- -teta 
}   
   
########################################################################################   
########################################################################################



if(BivD=="N"){


c.copula.be1 <- pnorm( (qnorm(p2) - teta*qnorm(p1))/sqrt(1 - teta^2)   )  # BiCopHfunc(p1, p2, family=1, par=teta)$hfunc1                          
c.copula.be2 <- pnorm( (qnorm(p1) - teta*qnorm(p2))/sqrt(1 - teta^2)   )  # BiCopHfunc(p1, p2, family=1, par=teta)$hfunc2

c.copula.theta <- dbinorm(qnorm(p1),qnorm(p2), cov12=teta)*(1/cosh(teta.st)^2) 

c.copula2.be1 <- dnorm((qnorm(p2)-teta*qnorm(p1))/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qnorm(p1)^2/2) # BiCopHfuncDeriv(p2, p1, 1, par=teta, deriv="u2")    
c.copula2.be2 <- dnorm((qnorm(p1)-teta*qnorm(p2))/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qnorm(p2)^2/2) # BiCopHfuncDeriv(p1, p2, 1, par=teta, deriv="u2")

c.copula2.be1be2 <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) # BiCopPDF(p1, p2, 1, par=teta)

c.copula2.be1th <- -(dnorm((qnorm(p2) - tanh(teta.st) * qnorm(p1))/sqrt(1 - tanh(teta.st)^2)) * 
     (1/cosh(teta.st)^2 * qnorm(p1)/sqrt(1 - tanh(teta.st)^2) - (qnorm(p2) - 
         tanh(teta.st) * qnorm(p1)) * (0.5 * (2 * (1/cosh(teta.st)^2 * 
         tanh(teta.st)) * (1 - tanh(teta.st)^2)^-0.5))/sqrt(1 - 
         tanh(teta.st)^2)^2)) 

c.copula2.be2th <- -(dnorm((qnorm(p1) - tanh(teta.st) * qnorm(p2))/sqrt(1 - tanh(teta.st)^2)) * 
    (1/cosh(teta.st)^2 * qnorm(p2)/sqrt(1 - tanh(teta.st)^2) - (qnorm(p1) - 
        tanh(teta.st) * qnorm(p2)) * (0.5 * (2 * (1/cosh(teta.st)^2 * 
        tanh(teta.st)) * (1 - tanh(teta.st)^2)^-0.5))/sqrt(1 - 
        tanh(teta.st)^2)^2))
 
bit1.th2 <- (2 * pi * (0.5 * (2 * (1/cosh(teta.st)^2 * tanh(teta.st)) * (1 - 
    tanh(teta.st)^2)^-0.5))/(2 * pi * sqrt(1 - tanh(teta.st)^2))^2 * 
    exp(-1/(2 * (1 - tanh(teta.st)^2)) * (qnorm(p1)^2 + qnorm(p2)^2 - 2 * 
        tanh(teta.st) * qnorm(p1) * qnorm(p2))) - 1/(2 * pi * sqrt(1 - 
    tanh(teta.st)^2)) * (exp(-1/(2 * (1 - tanh(teta.st)^2)) * 
    (qnorm(p1)^2 + qnorm(p2)^2 - 2 * tanh(teta.st) * qnorm(p1) * qnorm(p2))) * (-1/(2 * 
    (1 - tanh(teta.st)^2)) * (2 * (1/cosh(teta.st)^2) * qnorm(p1) * 
    qnorm(p2)) + 2 * (2 * (1/cosh(teta.st)^2 * tanh(teta.st)))/(2 * 
    (1 - tanh(teta.st)^2))^2 * (qnorm(p1)^2 + qnorm(p2)^2 - 2 * tanh(teta.st) * 
    qnorm(p1) * qnorm(p2)))))/cosh(teta.st)^2 - 1/(2 * pi * sqrt(1 - tanh(teta.st)^2)) * 
    exp(-1/(2 * (1 - tanh(teta.st)^2)) * (qnorm(p1)^2 + qnorm(p2)^2 - 2 * 
        tanh(teta.st) * qnorm(p1) * qnorm(p2))) * 1 * (2 * (sinh(teta.st) * 
    cosh(teta.st)))/(cosh(teta.st)^2)^2 


bit1.th2ATE <- 0.5 * (pi * (0.5 * (2 * teta * (1 - teta^2)^-0.5)))/(pi * sqrt(1 - 
    teta^2))^2 * (exp(-0.5/(1 - teta^2) * (qnorm(p1)^2 + qnorm(p2)^2 - 2 * 
    teta * qnorm(p1) * qnorm(p2)))) - 0.5/(pi * sqrt(1 - teta^2)) * (exp(-0.5/(1 - 
    teta^2) * (qnorm(p1)^2 + qnorm(p2)^2 - 2 * teta * qnorm(p1) * qnorm(p2))) * (-0.5/(1 - 
    teta^2) * (2 * qnorm(p1) * qnorm(p2)) + 0.5 * (2 * teta)/(1 - teta^2)^2 * 
    (qnorm(p1)^2 + qnorm(p2)^2 - 2 * teta * qnorm(p1) * qnorm(p2))))




}







if(BivD=="F"){

# 1 - exp(-teta) = -expm1(-teta)
# recall log1p as well if needed. 
# epcl <- -expm1(-teta) 


  c.copula.be1 <-     (exp(teta)* (-1 + exp(p2* teta)))/(-exp((p1 + p2)* teta) + 
                       exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))


 
  c.copula.be2 <-   (exp(teta)* (-1 + exp(p1* teta)))/(-exp((p1 + p2)* teta) + 
                     exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))


  c.copula.theta <-   (exp(teta)* (1/(-1 + exp(teta)) + (-1 - exp(p2* teta)* (-1 + p1) + p1 - 
     exp(p1* teta)* (-1 + p2) + p2)/(
    exp((p1 + p2)* teta) - exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta))))* teta + 
 log((exp(-(p1 + p2)* teta)* (-exp((p1 + p2)* teta) + 
     exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta))))/(-1 + exp(teta))))/teta^2

       
c.copula2.be1 <-   (exp(teta + 
  p1* teta)* (-1 + exp(p2* teta))* (-exp(teta) + exp(
   p2* teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2
  
 
 c.copula2.be2 <-     (exp(teta + 
  p2* teta)* (-1 + exp(p1* teta))* (-exp(teta) + exp(
   p1* teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))^2


c.copula2.be1be2 <- (exp((1 + p1 + p2)* teta)* (-1 + exp(teta))* teta)/(exp((p1 + p2)* teta) - 
  exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2


c.copula2.be1th <- (exp(teta + 
  p1 *teta)* (exp(2* p2* teta)* (-1 + p1) + exp(teta)* p1 - 
   exp(p2* teta)* (-1 + p1 + exp(teta)* (p1 - p2) + p2)))/(exp((p1 + 
     p2)* teta) - exp(teta)* (-1 + exp(p1* teta) + exp(p2* teta)))^2

c.copula2.be2th <-  (exp(teta + 
  p2 *teta)* (exp(2* p1* teta)* (-1 + p2) + exp(teta)* p2 - 
   exp(p1* teta)* (-1 + p2 + exp(teta)* (p2 - p1) + p1)))/(exp((p2 + 
     p1)* teta) - exp(teta)* (-1 + exp(p2* teta) + exp(p1* teta)))^2

#####
bit1.th2 <- bit1.th2ATE <- 1/teta^2 * ((1/(1 - exp(-teta)) * (exp(-teta) - (exp(-teta * 
    p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))) - 2 * teta/(teta^2)^2 * 
    log(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * 
        p1)) * (1 - exp(-teta * p2)))) + (1/teta^2 * ((1/(1 - 
    exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - exp(-teta)) * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))) - -1/teta * ((1/(1 - exp(-teta)) * (exp(-teta) + 
    (exp(-teta * p1) * p1 * (exp(-teta * p2) * p2) - exp(-teta * 
        p1) * p1 * p1 * (1 - exp(-teta * p2)) + (exp(-teta * 
        p1) * p1 * (exp(-teta * p2) * p2) - (1 - exp(-teta * 
        p1)) * (exp(-teta * p2) * p2 * p2)))) + exp(-teta)/(1 - 
    exp(-teta))^2 * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) + (exp(-teta)/(1 - exp(-teta))^2 * (exp(-teta) - (exp(-teta * 
    p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * p2))) - (exp(-teta)/(1 - exp(-teta))^2 + 
    exp(-teta) * (2 * (exp(-teta) * (1 - exp(-teta))))/((1 - 
        exp(-teta))^2)^2) * ((1 - exp(-teta)) - (1 - exp(-teta * 
    p1)) * (1 - exp(-teta * p2)))))/(1/(1 - exp(-teta)) * ((1 - 
    exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) + 
    (1/(1 - exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * 
        (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * 
        p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) * (1/(1 - 
        exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
        exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * 
        p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - 
        exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * 
        (1 - exp(-teta * p2))))^2))


}





if(BivD %in% c("C0","C90","C180","C270")){



c.copula.be1 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p1^((-teta) - 1) * (-teta)))

c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))
  
  
  c.copula.thet <- ((-1 + p1^-teta + p2^-teta)^(-1/
  teta) *((teta *(p2^teta *log(p1) + p1^teta* log(p2)))/(
   p2^teta - p1^teta* (-1 + p2^teta)) + 
   log(-1 + p1^-teta + p2^-teta)))/teta^2
  
  c.copula.theta <- c.copula.thet*derteta.derteta.st


c.copula2.be1 <- (p1^(-2 + teta)* p2^teta* (-1 + p1^-teta + p2^-teta)^(-1/
  teta)* (-1 + p2^teta)* (1 + teta))/(p2^teta - 
  p1^teta* (-1 + p2^teta))^2
 

 c.copula2.be2 <- (p2^(-2 + teta)* p1^teta* (-1 + p2^-teta + p1^-teta)^(-1/
  teta)* (-1 + p1^teta)* (1 + teta))/(p1^teta - 
  p2^teta* (-1 + p1^teta))^2


c.copula2.be1be2 <- p1^(-1 - teta)* p2^(-1 - teta)* (-1 + p1^-teta + p2^-teta)^(-2 - 1/
  teta) *(1 + teta)
 
  
c.copula2.be1th <- (p2^teta *(-1 + p1^-teta + p2^-teta)^(-1/
  teta)* (teta* (p2^teta + p1^teta* (-1 + p2^teta)* teta)* log(p1) + 
   p1^teta* teta* (1 + teta)* log(
     p2) - (-p1^teta + (-1 + p1^teta)* p2^teta)* log(-1 + p1^-teta + 
      p2^-teta)))/(p1* (p2^teta - p1^teta* (-1 + p2^teta))^2* teta)


c.copula2.be2th <- (p1^teta *(-1 + p2^-teta + p1^-teta)^(-1/
  teta)* (teta* (p1^teta + p2^teta* (-1 + p1^teta)* teta)* log(p2) + 
   p2^teta* teta* (1 + teta)* log(
     p1) - (-p2^teta + (-1 + p2^teta)* p1^teta)* log(-1 + p2^-teta + 
      p1^-teta)))/(p2* (p1^teta - p2^teta* (-1 + p1^teta))^2* teta)
    
 
bit1.th2ATE <-(1/((p2^teta - p1^teta* (-1 + p2^teta))^2* teta^4))*(-1 + p1^-teta + 
      p2^-teta)^(-1/teta)* (p2^
       teta *teta^2* (p2^teta + p1^teta* (-1 + p2^teta)* teta)* log(p1)^2 + 
      p1^teta* teta^2* (p1^teta + (-1 + p1^teta)* p2^teta* teta)* log(p2)^2 + 
      2 *p2^teta* teta* log(p1)* (p1^teta* teta* (1 + teta)* log(
           p2) + (-p1^teta + (-1 + p1^teta)* p2^teta)* (teta - 
            log(-1 + p1^-teta + p2^-teta))) + 
      2* p1^teta* (-p1^teta + (-1 + p1^teta)* p2^teta)* teta* log(
        p2) *(teta - log(-1 + p1^-teta + p2^-teta)) - (p1^teta + p2^teta -
          p1^teta* p2^teta)^2* (2* teta - 
      log(-1 + p1^-teta + p2^-teta))* log(-1 + p1^-teta + p2^-teta))  
   
    

bit1.th2 <- bit1.th2ATE*derteta.derteta.st^2 + c.copula.thet*der2teta.derteta.stteta.st  




}








if(BivD %in% c("G0","G90","G180","G270")){



  c.copula.be1 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta)* ((-log(p1))^teta + (-log(p2))^teta)^(-1 + 1/teta))/p1


  c.copula.be2 <- (exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/teta)))* (-log(p2))^(-1 + 
  teta)* ((-log(p2))^teta + (-log(p1))^teta)^(-1 + 1/teta))/p2


  c.copula.thet <-   (1/(teta^2))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* ((-log(p1))^teta + (-log(p2))^teta)^(-1 + 1/
  teta)* (-teta* (-log(p1))^
    teta* log(-log(p1)) + ((-log(p1))^teta + (-log(p2))^
      teta)* log((-log(p1))^teta + (-log(p2))^teta) - 
   teta *(-log(p2))^teta* log(-log(p2)))

  c.copula.theta <- c.copula.thet*derteta.derteta.st
    
    
    
   
c.copula2.be1 <-(1/(p1^2))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* (-log(p1))^(-2 + 
  teta) *((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* ((-log(p1))^
    teta* (log(p1) + ((-log(p1))^teta + (-log(p2))^teta)^(1/
      teta)) + (1 - teta + log(p1))* (-log(p2))^teta)


                  
c.copula2.be2 <- (1/(p2^2))*exp(-((-log(p2))^teta + (-log(p1))^teta)^((1/
  teta)))* (-log(p2))^(-2 + 
  teta) *((-log(p2))^teta + (-log(p1))^teta)^(-2 + 1/
  teta)* ((-log(p2))^
    teta* (log(p2) + ((-log(p2))^teta + (-log(p1))^teta)^(1/
      teta)) + (1 - teta + log(p2))* (-log(p1))^teta)


c.copula2.be1be2 <- (exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/teta)))* (-log(p1))^(-1 + 
  teta) *(-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(1/
   teta))* ((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* (-log(p2))^(-1 + teta))/(p1 *p2)



c.copula2.be1th <-(1/(p1 *teta^2 *log(p1)))*exp(
 teta.st - ((-log(p1))^teta + (-log(p2))^teta)^(1/
  teta))* (-log(p1))^teta* ((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* (-teta* (-(-log(p1))^
        teta *(-1 + ((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) + 
      teta* (-log(p2))^teta) *log(-log(
       p1)) - (-1 + ((-log(p1))^teta + (-log(p2))^teta)^(1/
      teta))* ((-log(p1))^teta + (-log(p2))^teta)* log((-log(p1))^
      teta + (-log(p2))^teta) + 
   teta* (-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(1/
      teta))* (-log(p2))^teta *log(-log(p2)))


c.copula2.be2th <-(1/(p2 *teta^2 *log(p2)))*exp(
 teta.st - ((-log(p2))^teta + (-log(p1))^teta)^(1/
  teta))* (-log(p2))^teta* ((-log(p2))^teta + (-log(p1))^teta)^(-2 + 1/
  teta)* (-teta* (-(-log(p2))^
        teta *(-1 + ((-log(p2))^teta + (-log(p1))^teta)^(1/teta)) + 
      teta* (-log(p1))^teta) *log(-log(
       p2)) - (-1 + ((-log(p2))^teta + (-log(p1))^teta)^(1/
      teta))* ((-log(p2))^teta + (-log(p1))^teta)* log((-log(p2))^
      teta + (-log(p1))^teta) + 
   teta* (-1 + teta + ((-log(p2))^teta + (-log(p1))^teta)^(1/
      teta))* (-log(p1))^teta *log(-log(p1)))



bit1.th2ATE <- (1/(teta^4))*exp(-((-log(p1))^teta + (-log(p2))^teta)^((1/
  teta)))* ((-log(p1))^teta + (-log(p2))^teta)^(-2 + 1/
  teta)* (teta^2* (-log(p1))^
    teta* ((-log(p1))^
       teta* (-1 + ((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) - 
      teta* (-log(p2))^teta) *log(-log(
       p1))^2 + (-1 + ((-log(p1))^teta + (-log(p2))^teta)^(1/
      teta))* ((-log(p1))^teta + (-log(p2))^teta)^2* log((-log(p1))^
      teta + (-log(p2))^teta)^2 - 
   2* teta* ((-log(p1))^teta + (-log(p2))^teta)* log((-log(p1))^
      teta + (-log(p2))^teta)* ((-log(p1))^
      teta + (-log(p2))^
       teta* (1 + (-1 + ((-log(p1))^teta + (-log(p2))^teta)^(1/
            teta))* log(-log(p2)))) + 
   teta^2* (-log(p2))^
    teta* log(-log(p2))* ((-log(p1))^
       teta* (2 - teta* log(-log(p2))) + (-log(p2))^
       teta* (2 + (-1 + ((-log(p1))^teta + (-log(p2))^teta)^(1/
            teta))* log(-log(p2)))) + 
   2* teta* (-log(p1))^
    teta* log(-log(p1))* (-(-1 + ((-log(p1))^teta + (-log(p2))^teta)^(1/
          teta))* ((-log(p1))^teta + (-log(p2))^teta)* log((-log(p1))^
         teta + (-log(p2))^teta) + 
      teta* ((-log(p1))^
         teta + (-log(p2))^
          teta *(1 + (-1 + teta + ((-log(p1))^teta + (-log(p2))^teta)^(
               1/teta))* log(-log(p2))))))

 
bit1.th2 <- bit1.th2ATE*derteta.derteta.st^2 + c.copula.thet*der2teta.derteta.stteta.st   

}






if(BivD %in% c("J0","J90","J180","J270")){


  c.copula.be1 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p1)^(teta - 1) * teta - (1 - p1)^(teta - 
    1) * teta * (1 - p2)^teta))


 
  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))




  c.copula.thet <-    (((1 - p1)^
   teta - (-1 + (1 - p1)^teta) *(1 - p2)^
    teta)^(1/teta)* (log((1 - p1)^
     teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta) + (
   teta* ((1 - p1)^
       teta* (-1 + (1 - p2)^teta)* log(
        1 - p1) + (-1 + (1 - p1)^teta) *(1 - p2)^
       teta* log(1 - p2)))/((1 - p1)^
    teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)))/teta^2
        
  c.copula.theta <- c.copula.thet*derteta.derteta.st        
  



c.copula2.be1 <- (1 - p1)^(-2 + 
  teta)* (-1 + (1 - p2)^teta)* ((1 - p1)^
   teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^teta* (-1 + teta)


               
c.copula2.be2 <- (1 - p2)^(-2 + 
  teta)* (-1 + (1 - p1)^teta)* ((1 - p2)^
   teta - (-1 + (1 - p2)^teta)* (1 - p1)^teta)^(-2 + 1/
  teta) *(1 - p1)^teta* (-1 + teta)



c.copula2.be1be2 <- (1 - p1)^(-1 + 
  teta)* ((1 - p1)^teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(1 - p2)^(-1 + 
  teta)* (-(-1 + (1 - p1)^teta)* (-1 + (1 - p2)^teta) + teta)


c.copula2.be1th <--(1/(teta^2))*
 exp(teta.st)* (1 - p1)^(-1 + 
   teta)* ((1 - p1)^teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
   teta)* ((-1 + (1 - p2)^teta)* teta *((1 - p1)^
       teta + (1 - p2)^teta *(-(1 - p1)^teta + teta))* log(
      1 - p1) + (-1 + (1 - p2)^
       teta)* (-(1 - p1)^teta + (-1 + (1 - p1)^teta)* (1 - p2)^
        teta)* log((1 - p1)^
       teta - (-1 + (1 - p1)^teta) *(1 - p2)^teta) - (1 - p2)^
     teta *((-1 + (1 - p1)^teta)* (-1 + (1 - p2)^teta) - teta)* teta* log(
      1 - p2))


c.copula2.be2th <--(1/(teta^2))*
 exp(teta.st)* (1 - p2)^(-1 + 
   teta)* ((1 - p2)^teta - (-1 + (1 - p2)^teta)* (1 - p1)^teta)^(-2 + 1/
   teta)* ((-1 + (1 - p1)^teta)* teta *((1 - p2)^
       teta + (1 - p1)^teta *(-(1 - p2)^teta + teta))* log(
      1 - p2) + (-1 + (1 - p1)^
       teta)* (-(1 - p2)^teta + (-1 + (1 - p2)^teta)* (1 - p1)^
        teta)* log((1 - p2)^
       teta - (-1 + (1 - p2)^teta) *(1 - p1)^teta) - (1 - p1)^
     teta *((-1 + (1 - p2)^teta)* (-1 + (1 - p1)^teta) - teta)* teta* log(
      1 - p1))
      
bit1.th2ATE <- (1/(teta^4))*((1 - p1)^
   teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^(-2 + 1/
  teta) *(-(1 - p1)^
     teta* (-1 + (1 - p2)^
      teta) *teta^2* ((1 - p1)^teta* (-1 + (1 - p2)^teta) - (1 - p2)^
       teta* teta)* log(
     1 - p1)^2 - ((1 - p1)^
      teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^2* log((1 - p1)^
      teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)^2 + 
   2 *(-(1 - p1)^teta + (-1 + (1 - p1)^teta)* (1 - p2)^
       teta)* teta* log((1 - p1)^
      teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta)* ((1 - p1)^
      teta + (-1 + (1 - p1)^teta)* (1 - p2)^
       teta* (-1 + log(1 - p2))) + (-1 + (1 - p1)^teta)* (1 - p2)^
    teta *teta^2 *log(
     1 - p2)* (-(-1 + (1 - p1)^teta)* (1 - p2)^
       teta* (-2 + log(1 - p2)) + (1 - p1)^
       teta* (-2 + teta* log(1 - p2))) + 
   2 *(1 - p1)^
    teta* teta* log(
     1 - p1)* ((-1 + (1 - p2)^
         teta)* (-(1 - p1)^teta + (-1 + (1 - p1)^teta)* (1 - p2)^
          teta)* log((1 - p1)^
         teta - (-1 + (1 - p1)^teta)* (1 - p2)^teta) + 
      teta* ((1 - p1)^
         teta - (-1 + (1 - p1)^teta)* (1 - p2)^(
          2* teta) *(-1 + log(1 - p2)) + (1 - p2)^
          teta *(1 - 2* (1 - p1)^
             teta + (-1 + (1 - p1)^teta + teta)* log(1 - p2)))))
 
 
bit1.th2 <- bit1.th2ATE*derteta.derteta.st^2 + c.copula.thet*der2teta.derteta.stteta.st   

}



if(BivD %in% c("C90","J90","G90") ) {

c.copula.be1     <- c.copula.be1
c.copula.be2     <- 1 - c.copula.be2
c.copula.theta   <- - c.copula.theta
c.copula2.be1    <- - c.copula2.be1
c.copula2.be2    <- - c.copula2.be2
c.copula2.be1be2 <- c.copula2.be1be2
c.copula2.be1th  <- c.copula2.be1th 
c.copula2.be2th  <- - c.copula2.be2th
bit1.th2ATE      <- - bit1.th2ATE  
bit1.th2         <- - bit1.th2 

}  

if(BivD %in% c("C180","J180","G180") ) {

c.copula.be1     <- 1 - c.copula.be1 
c.copula.be2     <- 1 - c.copula.be2
c.copula.theta   <- c.copula.theta
c.copula2.be1    <- c.copula2.be1 
c.copula2.be2    <- c.copula2.be2
c.copula2.be1be2 <- c.copula2.be1be2
c.copula2.be1th  <- - c.copula2.be1th
c.copula2.be2th  <- - c.copula2.be2th
bit1.th2ATE      <- bit1.th2ATE
bit1.th2         <- bit1.th2  


}  


if(BivD %in% c("C270","J270","G270") ) {

c.copula.be1     <- 1 - c.copula.be1
c.copula.be2     <- c.copula.be2
c.copula.theta   <- - c.copula.theta
c.copula2.be1    <- - c.copula2.be1
c.copula2.be2    <- - c.copula2.be2
c.copula2.be1be2 <- c.copula2.be1be2
c.copula2.be1th  <- - c.copula2.be1th
c.copula2.be2th  <-   c.copula2.be2th
bit1.th2ATE      <- - bit1.th2ATE
bit1.th2         <- - bit1.th2

}   




epsilon <- 0.0000001 
max.p   <- 0.9999999
  
c.copula.be2 <- ifelse(c.copula.be2 > max.p, max.p, c.copula.be2) 
c.copula.be2 <- ifelse(c.copula.be2 < epsilon,     epsilon, c.copula.be2)
c.copula.be1 <- ifelse(c.copula.be1 > max.p, max.p, c.copula.be1) 
c.copula.be1 <- ifelse(c.copula.be1 < epsilon,     epsilon, c.copula.be1)





list(
c.copula.be1     = c.copula.be1,    
c.copula.be2     = c.copula.be2,    
c.copula.theta   = c.copula.theta,  
c.copula2.be1    = c.copula2.be1,  
c.copula2.be2    = c.copula2.be2,   
c.copula2.be1be2 = c.copula2.be1be2,
c.copula2.be1th  = c.copula2.be1th, 
c.copula2.be2th  = c.copula2.be2th, 
bit1.th2ATE      = bit1.th2ATE,     
bit1.th2         = bit1.th2 )     


}


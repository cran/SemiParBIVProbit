copgHs <- function(p1,p2,teta,teta.st,BivD,nC,nu){

epsilon <- .Machine$double.eps*10^6

if(BivD=="N"){

#c.copula.be1 <- pnorm( (eta2-teta*eta1)*(1/sqrt( pmax(10000*.Machine$double.eps, 1-teta^2) )) )
#c.copula.be2 <- pnorm( (eta1-teta*eta2)*(1/sqrt( pmax(10000*.Machine$double.eps, 1-teta^2) )) )
#c.copula.theta <- dbinorm(eta1,eta2, cov12=teta)*(1/cosh(teta.st)^2) 


c.copula.be1 <- BiCopHfunc(p1, p2, family=1, par=teta)$hfunc1
c.copula.be2 <- BiCopHfunc(p1, p2, family=1, par=teta)$hfunc2
c.copula.theta <- dbinorm(qnorm(p1),qnorm(p2), cov12=teta)*(1/cosh(teta.st)^2) 




#c.copula2.be1 <- dnorm((eta2-teta*eta1)*(1/sqrt( pmax(10000*.Machine$double.eps, 1-teta^2) )))*-teta*(1/sqrt( pmax(10000*.Machine$double.eps, 1-teta^2) ))*sqrt(2*pi)/exp(-eta1^2/2)  
#c.copula2.be2 <- dnorm((eta1-teta*eta2)*(1/sqrt( pmax(10000*.Machine$double.eps, 1-teta^2) )))*-teta*(1/sqrt( pmax(10000*.Machine$double.eps, 1-teta^2) ))*sqrt(2*pi)/exp(-eta2^2/2)



c.copula2.be1 <- BiCopHfuncDeriv(p2, p1, 1, par=teta, deriv="u2")                
c.copula2.be2 <- BiCopHfuncDeriv(p1, p2, 1, par=teta, deriv="u2")


#c.copula2.be1be2 <- (1/sqrt( pmax(10000*.Machine$double.eps, 1-teta^2) ))*exp( -(  teta^2*(eta1^2+eta2^2)-2*teta*eta1*eta2     )/(2*(1-teta^2))  )


c.copula2.be1be2 <- BiCopPDF(p1, p2, 1, par=teta)




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


}



if(BivD=="T"){

c.copula.be1 <- BiCopHfunc(p1, p2, family=2, par=teta, par2=nu)$hfunc1
c.copula.be2 <- BiCopHfunc(p1, p2, family=2, par=teta, par2=nu)$hfunc2
c.copula.theta <- ( 1 + (qt(p1,nu)^2+qt(p2,nu)^2-2*teta*qt(p1,nu)*qt(p2,nu))/(nu*(1-teta^2)) )^(-nu/2)/(2*pi*sqrt(1-teta^2))*(1/cosh(teta.st)^2) 

c.copula2.be1    <- BiCopHfuncDeriv(p2, p1, 2, par=teta, par2=nu, deriv="u2")  
c.copula2.be2    <- BiCopHfuncDeriv(p1, p2, 2, par=teta, par2=nu, deriv="u2")
c.copula2.be1be2 <- BiCopPDF(p1, p2, 2, par=teta, par2=nu)
c.copula2.be1th  <- BiCopHfuncDeriv(p2, p1, 2, par=teta, par2=nu, deriv="par")/cosh(teta.st)^2 
c.copula2.be2th  <- BiCopHfuncDeriv(p1, p2, 2, par=teta, par2=nu, deriv="par")/cosh(teta.st)^2
 

bit1.th2 <- -(((1 + (qt(p1,nu)^2 + qt(p2,nu)^2 - 2 * tanh(teta.st) * qt(p1,nu) * qt(p2,nu))/(nu * 
    (1 - tanh(teta.st)^2)))^((-nu/2) - 1) * ((-nu/2) * (2 * (1/cosh(teta.st)^2) * 
    qt(p1,nu) * qt(p2,nu)/(nu * (1 - tanh(teta.st)^2)) - (qt(p1,nu)^2 + qt(p2,nu)^2 - 
    2 * tanh(teta.st) * qt(p1,nu) * qt(p2,nu)) * (nu * (2 * (1/cosh(teta.st)^2 * 
    tanh(teta.st))))/(nu * (1 - tanh(teta.st)^2))^2))/(2 * pi * 
    sqrt(1 - tanh(teta.st)^2)) - (1 + (qt(p1,nu)^2 + qt(p2,nu)^2 - 2 * 
    tanh(teta.st) * qt(p1,nu) * qt(p2,nu))/(nu * (1 - tanh(teta.st)^2)))^(-nu/2) * 
    (2 * pi * (0.5 * (2 * (1/cosh(teta.st)^2 * tanh(teta.st)) * 
        (1 - tanh(teta.st)^2)^-0.5)))/(2 * pi * sqrt(1 - tanh(teta.st)^2))^2)/cosh(teta.st)^2 + 
    (1 + (qt(p1,nu)^2 + qt(p2,nu)^2 - 2 * tanh(teta.st) * qt(p1,nu) * qt(p2,nu))/(nu * 
        (1 - tanh(teta.st)^2)))^(-nu/2)/(2 * pi * sqrt(1 - tanh(teta.st)^2)) * 
        1 * (2 * (sinh(teta.st) * cosh(teta.st)))/(cosh(teta.st)^2)^2) 


}






if(BivD=="C0"){

c.copula.be1 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p1^((-teta) - 1) * (-teta)))
  c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))
  c.copula.theta <- ((p1^(-teta) + p2^(-teta) - 1)^(-1/teta) * (log((p1^(-teta) + 
    p2^(-teta) - 1)) * (1/teta^2)) - (p1^(-teta) + p2^(-teta) - 
    1)^((-1/teta) - 1) * ((-1/teta) * (p2^(-teta) * log(p2) + 
    p1^(-teta) * log(p1))))*exp(teta.st)

c.copula2.be1 <-  (p1^(-teta) + p2^(-teta) - 1)^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (p1^((-teta) - 1) * (-teta))) * ((-1/teta) * (p1^((-teta) - 
    1) * (-teta))) + (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 
    1) * ((-1/teta) * (p1^(((-teta) - 1) - 1) * ((-teta) - 1) * 
    (-teta)))
                 
 c.copula2.be2 <- (p1^(-teta) + p2^(-teta) - 1)^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (p2^((-teta) - 1) * (-teta))) * ((-1/teta) * (p2^((-teta) - 
    1) * (-teta))) + (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 
    1) * ((-1/teta) * (p2^(((-teta) - 1) - 1) * ((-teta) - 1) * 
    (-teta)))



c.copula2.be1be2 <- (p1^(-teta) + p2^(-teta) - 1)^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (p2^((-teta) - 1) * (-teta))) * ((-1/teta) * (p1^((-teta) - 
    1) * (-teta)))


c.copula2.be1th <- ((p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 
    1)^((-1/(exp(teta.st) + epsilon)) - 1) * (log((p1^(-(exp(teta.st) + 
    epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2)) - (p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^(((-1/(exp(teta.st) + epsilon)) - 1) - 1) * 
    (((-1/(exp(teta.st) + epsilon)) - 1) * (p2^(-(exp(teta.st) + 
        epsilon)) * (log(p2) * exp(teta.st)) + p1^(-(exp(teta.st) + 
        epsilon)) * (log(p1) * exp(teta.st))))) * ((-1/(exp(teta.st) + 
    epsilon)) * (p1^((-(exp(teta.st) + epsilon)) - 1) * (-(exp(teta.st) + 
    epsilon)))) + (p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^((-1/(exp(teta.st) + epsilon)) - 1) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2 * (p1^((-(exp(teta.st) + epsilon)) - 1) * (-(exp(teta.st) + 
    epsilon))) - (-1/(exp(teta.st) + epsilon)) * (p1^((-(exp(teta.st) + 
    epsilon)) - 1) * exp(teta.st) + p1^((-(exp(teta.st) + epsilon)) - 
    1) * (log(p1) * exp(teta.st)) * (-(exp(teta.st) + epsilon))))


c.copula2.be2th <- ((p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 
    1)^((-1/(exp(teta.st) + epsilon)) - 1) * (log((p1^(-(exp(teta.st) + 
    epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2)) - (p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^(((-1/(exp(teta.st) + epsilon)) - 1) - 1) * 
    (((-1/(exp(teta.st) + epsilon)) - 1) * (p2^(-(exp(teta.st) + 
        epsilon)) * (log(p2) * exp(teta.st)) + p1^(-(exp(teta.st) + 
        epsilon)) * (log(p1) * exp(teta.st))))) * ((-1/(exp(teta.st) + 
    epsilon)) * (p2^((-(exp(teta.st) + epsilon)) - 1) * (-(exp(teta.st) + 
    epsilon)))) + (p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^((-1/(exp(teta.st) + epsilon)) - 1) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2 * (p2^((-(exp(teta.st) + epsilon)) - 1) * (-(exp(teta.st) + 
    epsilon))) - (-1/(exp(teta.st) + epsilon)) * (p2^((-(exp(teta.st) + 
    epsilon)) - 1) * exp(teta.st) + p2^((-(exp(teta.st) + epsilon)) - 
    1) * (log(p2) * exp(teta.st)) * (-(exp(teta.st) + epsilon))))



bit1.th2 <-((p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 
    1)^(-1/(exp(teta.st) + epsilon)) * (log((p1^(-(exp(teta.st) + 
    epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2)) - (p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^((-1/(exp(teta.st) + epsilon)) - 1) * ((-1/(exp(teta.st) + 
    epsilon)) * (p2^(-(exp(teta.st) + epsilon)) * (log(p2) * 
    exp(teta.st)) + p1^(-(exp(teta.st) + epsilon)) * (log(p1) * 
    exp(teta.st))))) * (log((p1^(-(exp(teta.st) + epsilon)) + 
    p2^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2)) + (p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^(-1/(exp(teta.st) + epsilon)) * (log((p1^(-(exp(teta.st) + 
    epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
    epsilon)))/((exp(teta.st) + epsilon)^2)^2) - (p2^(-(exp(teta.st) + 
    epsilon)) * (log(p2) * exp(teta.st)) + p1^(-(exp(teta.st) + 
    epsilon)) * (log(p1) * exp(teta.st)))/(p1^(-(exp(teta.st) + 
    epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 1) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2)) - (((p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^((-1/(exp(teta.st) + epsilon)) - 1) * (log((p1^(-(exp(teta.st) + 
    epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2)) - (p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^(((-1/(exp(teta.st) + epsilon)) - 1) - 1) * 
    (((-1/(exp(teta.st) + epsilon)) - 1) * (p2^(-(exp(teta.st) + 
        epsilon)) * (log(p2) * exp(teta.st)) + p1^(-(exp(teta.st) + 
        epsilon)) * (log(p1) * exp(teta.st))))) * ((-1/(exp(teta.st) + 
    epsilon)) * (p2^(-(exp(teta.st) + epsilon)) * (log(p2) * 
    exp(teta.st)) + p1^(-(exp(teta.st) + epsilon)) * (log(p1) * 
    exp(teta.st)))) + (p1^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
    epsilon)) - 1)^((-1/(exp(teta.st) + epsilon)) - 1) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2 * (p2^(-(exp(teta.st) + epsilon)) * (log(p2) * 
    exp(teta.st)) + p1^(-(exp(teta.st) + epsilon)) * (log(p1) * 
    exp(teta.st))) + (-1/(exp(teta.st) + epsilon)) * (p2^(-(exp(teta.st) + 
    epsilon)) * (log(p2) * exp(teta.st)) - p2^(-(exp(teta.st) + 
    epsilon)) * (log(p2) * exp(teta.st)) * (log(p2) * exp(teta.st)) + 
    (p1^(-(exp(teta.st) + epsilon)) * (log(p1) * exp(teta.st)) - 
        p1^(-(exp(teta.st) + epsilon)) * (log(p1) * exp(teta.st)) * 
            (log(p1) * exp(teta.st))))))


}
 
if(BivD=="C90"){

c.copula.be1 <- ((1 - p1)^(-(-teta)) + p2^(-(-teta)) - 1)^((-1/(-teta)) - 1) * 
    ((-1/(-teta)) * ((1 - p1)^((-(-teta)) - 1) * (-(-teta))))
  c.copula.be2 <- 1 - ((1 - p1)^(-(-teta)) + p2^(-(-teta)) - 1)^((-1/(-teta)) - 
    1) * ((-1/(-teta)) * (p2^((-(-teta)) - 1) * (-(-teta))))
  c.copula.theta <- (-(((1 - p1)^teta + p2^teta - 1)^((1/teta) - 1) * ((1/teta) * 
    ((1 - p1)^teta * log((1 - p1)) + p2^teta * log(p2))) - ((1 - 
    p1)^teta + p2^teta - 1)^(1/teta) * (log(((1 - p1)^teta + 
    p2^teta - 1)) * (1/teta^2))))*(-exp(teta.st))

  c.copula2.be1 <-  -(((1 - p1)^(-(-teta)) + p2^(-(-teta)) - 1)^((-1/(-teta)) - 1) * 
     ((-1/(-teta)) * ((1 - p1)^(((-(-teta)) - 1) - 1) * ((-(-teta)) - 
         1) * (-(-teta)))) + ((1 - p1)^(-(-teta)) + p2^(-(-teta)) - 
     1)^(((-1/(-teta)) - 1) - 1) * (((-1/(-teta)) - 1) * ((1 - 
     p1)^((-(-teta)) - 1) * (-(-teta)))) * ((-1/(-teta)) * ((1 - 
     p1)^((-(-teta)) - 1) * (-(-teta)))))
 
 
                   
  c.copula2.be2 <- -(((1 - p1)^(-(-teta)) + p2^(-(-teta)) - 1)^(((-1/(-teta)) - 
     1) - 1) * (((-1/(-teta)) - 1) * (p2^((-(-teta)) - 1) * (-(-teta)))) * 
     ((-1/(-teta)) * (p2^((-(-teta)) - 1) * (-(-teta)))) + ((1 - 
     p1)^(-(-teta)) + p2^(-(-teta)) - 1)^((-1/(-teta)) - 1) * 
     ((-1/(-teta)) * (p2^(((-(-teta)) - 1) - 1) * ((-(-teta)) - 
         1) * (-(-teta)))))
 
 
 c.copula2.be1be2 <- ((1 - p1)^(-(-teta)) + p2^(-(-teta)) - 1)^(((-1/(-teta)) - 1) - 
     1) * (((-1/(-teta)) - 1) * (p2^((-(-teta)) - 1) * (-(-teta)))) * 
     ((-1/(-teta)) * ((1 - p1)^((-(-teta)) - 1) * (-(-teta))))
 
 
 
 c.copula2.be1th <- ((((1 - p1)^teta + p2^teta - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
     1) * ((1 - p1)^teta * log((1 - p1)) + p2^teta * log(p2))) - 
     ((1 - p1)^teta + p2^teta - 1)^((1/teta) - 1) * (log(((1 - 
         p1)^teta + p2^teta - 1)) * (1/teta^2))) * ((1/teta) * 
     ((1 - p1)^(teta - 1) * teta)) + ((1 - p1)^teta + p2^teta - 
     1)^((1/teta) - 1) * ((1/teta) * ((1 - p1)^(teta - 1) * log((1 - 
     p1)) * teta + (1 - p1)^(teta - 1)) - 1/teta^2 * ((1 - p1)^(teta - 
     1) * teta)))*(-exp(teta.st))
   
 c.copula2.be2th <- (-((((1 - p1)^teta + p2^teta - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
     1) * ((1 - p1)^teta * log((1 - p1)) + p2^teta * log(p2))) - 
     ((1 - p1)^teta + p2^teta - 1)^((1/teta) - 1) * (log(((1 - 
         p1)^teta + p2^teta - 1)) * (1/teta^2))) * ((1/teta) * 
     (p2^(teta - 1) * teta)) + ((1 - p1)^teta + p2^teta - 1)^((1/teta) - 
     1) * ((1/teta) * (p2^(teta - 1) * log(p2) * teta + p2^(teta - 
     1)) - 1/teta^2 * (p2^(teta - 1) * teta))))*(-exp(teta.st))
 
 bit1.th2 <- -((((1 - p1)^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
     epsilon)) - 1)^(1/(-(exp(teta.st) + epsilon))) * (log(((1 - 
     p1)^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 
     1)) * (exp(teta.st)/(-(exp(teta.st) + epsilon))^2)) - ((1 - 
     p1)^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 
     1)^((1/(-(exp(teta.st) + epsilon))) - 1) * ((1/(-(exp(teta.st) + 
     epsilon))) * (p2^(-(exp(teta.st) + epsilon)) * (log(p2) * 
     exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + epsilon)) * (log((1 - 
     p1)) * exp(teta.st))))) * (log(((1 - p1)^(-(exp(teta.st) + 
     epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(-(exp(teta.st) + 
     epsilon))^2)) + ((1 - p1)^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
     epsilon)) - 1)^(1/(-(exp(teta.st) + epsilon))) * (log(((1 - 
     p1)^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + epsilon)) - 
     1)) * (exp(teta.st)/(-(exp(teta.st) + epsilon))^2 + exp(teta.st) * 
     (2 * (exp(teta.st) * (-(exp(teta.st) + epsilon))))/((-(exp(teta.st) + 
     epsilon))^2)^2) - (p2^(-(exp(teta.st) + epsilon)) * (log(p2) * 
     exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + epsilon)) * (log((1 - 
     p1)) * exp(teta.st)))/((1 - p1)^(-(exp(teta.st) + epsilon)) + 
     p2^(-(exp(teta.st) + epsilon)) - 1) * (exp(teta.st)/(-(exp(teta.st) + 
     epsilon))^2)) - ((((1 - p1)^(-(exp(teta.st) + epsilon)) + 
     p2^(-(exp(teta.st) + epsilon)) - 1)^((1/(-(exp(teta.st) + 
     epsilon))) - 1) * (log(((1 - p1)^(-(exp(teta.st) + epsilon)) + 
     p2^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(-(exp(teta.st) + 
     epsilon))^2)) - ((1 - p1)^(-(exp(teta.st) + epsilon)) + p2^(-(exp(teta.st) + 
     epsilon)) - 1)^(((1/(-(exp(teta.st) + epsilon))) - 1) - 1) * 
     (((1/(-(exp(teta.st) + epsilon))) - 1) * (p2^(-(exp(teta.st) + 
         epsilon)) * (log(p2) * exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + 
         epsilon)) * (log((1 - p1)) * exp(teta.st))))) * ((1/(-(exp(teta.st) + 
     epsilon))) * (p2^(-(exp(teta.st) + epsilon)) * (log(p2) * 
     exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + epsilon)) * (log((1 - 
     p1)) * exp(teta.st)))) + ((1 - p1)^(-(exp(teta.st) + epsilon)) + 
     p2^(-(exp(teta.st) + epsilon)) - 1)^((1/(-(exp(teta.st) + 
     epsilon))) - 1) * (exp(teta.st)/(-(exp(teta.st) + epsilon))^2 * 
     (p2^(-(exp(teta.st) + epsilon)) * (log(p2) * exp(teta.st)) + 
         (1 - p1)^(-(exp(teta.st) + epsilon)) * (log((1 - p1)) * 
             exp(teta.st))) + (1/(-(exp(teta.st) + epsilon))) * 
     (p2^(-(exp(teta.st) + epsilon)) * (log(p2) * exp(teta.st)) - 
         p2^(-(exp(teta.st) + epsilon)) * (log(p2) * exp(teta.st)) * 
             (log(p2) * exp(teta.st)) + ((1 - p1)^(-(exp(teta.st) + 
         epsilon)) * (log((1 - p1)) * exp(teta.st)) - (1 - p1)^(-(exp(teta.st) + 
         epsilon)) * (log((1 - p1)) * exp(teta.st)) * (log((1 - 
         p1)) * exp(teta.st)))))))
 


}



if(BivD=="C180"){

c.copula.be1 <- 1 - ((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^((-1/teta) - 1) * 
                  ((-1/teta) * ((1 - p1)^((-teta) - 1) * (-teta)))



  c.copula.be2 <- 1 - ((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^((-1/teta) - 1) * 
                  ((-1/teta) * ((1 - p2)^((-teta) - 1) * (-teta)))

  c.copula.theta <- (((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^(-1/teta) * (log(((1 - 
    p1)^(-teta) + (1 - p2)^(-teta) - 1)) * (1/teta^2)) - ((1 - 
    p1)^(-teta) + (1 - p2)^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * 
    ((1 - p2)^(-teta) * log((1 - p2)) + (1 - p1)^(-teta) * log((1 - 
        p1)))))*exp(teta.st)


 c.copula2.be1 <-  ((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * 
                   ((1 - p1)^(((-teta) - 1) - 1) * ((-teta) - 1) * (-teta))) + 
                   ((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^(((-1/teta) - 1) - 
                   1) * (((-1/teta) - 1) * ((1 - p1)^((-teta) - 1) * (-teta))) * 
                   ((-1/teta) * ((1 - p1)^((-teta) - 1) * (-teta)))
                  
 c.copula2.be2 <- ((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * 
                  ((1 - p2)^(((-teta) - 1) - 1) * ((-teta) - 1) * (-teta))) + 
                  ((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^(((-1/teta) - 1) - 
                  1) * (((-1/teta) - 1) * ((1 - p2)^((-teta) - 1) * (-teta))) * 
                  ((-1/teta) * ((1 - p2)^((-teta) - 1) * (-teta)))
                  

c.copula2.be1be2 <- ((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^(((-1/teta) - 1) - 
                    1) * (((-1/teta) - 1) * ((1 - p2)^((-teta) - 1) * (-teta))) * 
                    ((-1/teta) * ((1 - p1)^((-teta) - 1) * (-teta)))

c.copula2.be1th <- -((((1 - p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
                    epsilon)) - 1)^((-1/(exp(teta.st) + epsilon)) - 1) * (log(((1 - 
    p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
    epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + epsilon)^2)) - 
    ((1 - p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
        epsilon)) - 1)^(((-1/(exp(teta.st) + epsilon)) - 1) - 
        1) * (((-1/(exp(teta.st) + epsilon)) - 1) * ((1 - p2)^(-(exp(teta.st) + 
        epsilon)) * (log((1 - p2)) * exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + 
        epsilon)) * (log((1 - p1)) * exp(teta.st))))) * ((-1/(exp(teta.st) + 
    epsilon)) * ((1 - p1)^((-(exp(teta.st) + epsilon)) - 1) * 
    (-(exp(teta.st) + epsilon)))) + ((1 - p1)^(-(exp(teta.st) + 
    epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)^((-1/(exp(teta.st) + 
    epsilon)) - 1) * (exp(teta.st)/(exp(teta.st) + epsilon)^2 * 
    ((1 - p1)^((-(exp(teta.st) + epsilon)) - 1) * (-(exp(teta.st) + 
        epsilon))) - (-1/(exp(teta.st) + epsilon)) * ((1 - p1)^((-(exp(teta.st) + 
    epsilon)) - 1) * exp(teta.st) + (1 - p1)^((-(exp(teta.st) + 
    epsilon)) - 1) * (log((1 - p1)) * exp(teta.st)) * (-(exp(teta.st) + 
    epsilon)))))


c.copula2.be2th <- -((((1 - p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
    epsilon)) - 1)^((-1/(exp(teta.st) + epsilon)) - 1) * (log(((1 - 
    p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
    epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + epsilon)^2)) - 
    ((1 - p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
        epsilon)) - 1)^(((-1/(exp(teta.st) + epsilon)) - 1) - 
        1) * (((-1/(exp(teta.st) + epsilon)) - 1) * ((1 - p2)^(-(exp(teta.st) + 
        epsilon)) * (log((1 - p2)) * exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + 
        epsilon)) * (log((1 - p1)) * exp(teta.st))))) * ((-1/(exp(teta.st) + 
    epsilon)) * ((1 - p2)^((-(exp(teta.st) + epsilon)) - 1) * 
    (-(exp(teta.st) + epsilon)))) + ((1 - p1)^(-(exp(teta.st) + 
    epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)^((-1/(exp(teta.st) + 
    epsilon)) - 1) * (exp(teta.st)/(exp(teta.st) + epsilon)^2 * 
    ((1 - p2)^((-(exp(teta.st) + epsilon)) - 1) * (-(exp(teta.st) + 
        epsilon))) - (-1/(exp(teta.st) + epsilon)) * ((1 - p2)^((-(exp(teta.st) + 
    epsilon)) - 1) * exp(teta.st) + (1 - p2)^((-(exp(teta.st) + 
    epsilon)) - 1) * (log((1 - p2)) * exp(teta.st)) * (-(exp(teta.st) + 
    epsilon)))))

bit1.th2 <-(((1 - p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
    epsilon)) - 1)^(-1/(exp(teta.st) + epsilon)) * (log(((1 - 
    p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
    epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + epsilon)^2)) - 
    ((1 - p1)^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
        epsilon)) - 1)^((-1/(exp(teta.st) + epsilon)) - 1) * 
        ((-1/(exp(teta.st) + epsilon)) * ((1 - p2)^(-(exp(teta.st) + 
            epsilon)) * (log((1 - p2)) * exp(teta.st)) + (1 - 
            p1)^(-(exp(teta.st) + epsilon)) * (log((1 - p1)) * 
            exp(teta.st))))) * (log(((1 - p1)^(-(exp(teta.st) + 
    epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)) * 
    (exp(teta.st)/(exp(teta.st) + epsilon)^2)) + ((1 - p1)^(-(exp(teta.st) + 
    epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)^(-1/(exp(teta.st) + 
    epsilon)) * (log(((1 - p1)^(-(exp(teta.st) + epsilon)) + 
    (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
    epsilon)))/((exp(teta.st) + epsilon)^2)^2) - ((1 - p2)^(-(exp(teta.st) + 
    epsilon)) * (log((1 - p2)) * exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + 
    epsilon)) * (log((1 - p1)) * exp(teta.st)))/((1 - p1)^(-(exp(teta.st) + 
    epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2)) - ((((1 - p1)^(-(exp(teta.st) + epsilon)) + 
    (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)^((-1/(exp(teta.st) + 
    epsilon)) - 1) * (log(((1 - p1)^(-(exp(teta.st) + epsilon)) + 
    (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(exp(teta.st) + 
    epsilon)^2)) - ((1 - p1)^(-(exp(teta.st) + epsilon)) + (1 - 
    p2)^(-(exp(teta.st) + epsilon)) - 1)^(((-1/(exp(teta.st) + 
    epsilon)) - 1) - 1) * (((-1/(exp(teta.st) + epsilon)) - 1) * 
    ((1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * 
        exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + epsilon)) * 
        (log((1 - p1)) * exp(teta.st))))) * ((-1/(exp(teta.st) + 
    epsilon)) * ((1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - 
    p2)) * exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + epsilon)) * 
    (log((1 - p1)) * exp(teta.st)))) + ((1 - p1)^(-(exp(teta.st) + 
    epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)^((-1/(exp(teta.st) + 
    epsilon)) - 1) * (exp(teta.st)/(exp(teta.st) + epsilon)^2 * 
    ((1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * 
        exp(teta.st)) + (1 - p1)^(-(exp(teta.st) + epsilon)) * 
        (log((1 - p1)) * exp(teta.st))) + (-1/(exp(teta.st) + 
    epsilon)) * ((1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - 
    p2)) * exp(teta.st)) - (1 - p2)^(-(exp(teta.st) + epsilon)) * 
    (log((1 - p2)) * exp(teta.st)) * (log((1 - p2)) * exp(teta.st)) + 
    ((1 - p1)^(-(exp(teta.st) + epsilon)) * (log((1 - p1)) * 
        exp(teta.st)) - (1 - p1)^(-(exp(teta.st) + epsilon)) * 
        (log((1 - p1)) * exp(teta.st)) * (log((1 - p1)) * exp(teta.st))))))




}



if(BivD=="C270"){

c.copula.be1 <- 1 - (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * ((1/teta) * 
    (p1^((teta) - 1) * (teta)))

  c.copula.be2 <- (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * ((1/teta) * 
    ((1 - p2)^((teta) - 1) * (teta)))




  c.copula.theta <- (-((p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * ((1/teta) * 
    (p1^(teta) * log(p1) + (1 - p2)^(teta) * log((1 - p2)))) - 
    (p1^(teta) + (1 - p2)^(teta) - 1)^(1/teta) * (log((p1^(teta) + 
        (1 - p2)^(teta) - 1)) * (1/teta^2))))*(-exp(teta.st))

 c.copula2.be1 <- -((p1^(teta) + (1 - p2)^(teta) - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
     1) * (p1^((teta) - 1) * (teta))) * ((1/teta) * (p1^((teta) - 
     1) * (teta))) + (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 
     1) * ((1/teta) * (p1^(((teta) - 1) - 1) * ((teta) - 1) * 
     (teta))))
 
                   
  c.copula2.be2 <- -((p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * ((1/teta) * 
     ((1 - p2)^(((teta) - 1) - 1) * ((teta) - 1) * (teta))) + 
     (p1^(teta) + (1 - p2)^(teta) - 1)^(((1/teta) - 1) - 1) * 
         (((1/teta) - 1) * ((1 - p2)^((teta) - 1) * (teta))) * 
         ((1/teta) * ((1 - p2)^((teta) - 1) * (teta))))
 
 
 
 c.copula2.be1be2 <- (p1^(teta) + (1 - p2)^(teta) - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
     1) * ((1 - p2)^((teta) - 1) * (teta))) * ((1/teta) * (p1^((teta) - 
     1) * (teta)))
 
 
 c.copula2.be1th <- (-(((p1^(teta) + (1 - p2)^(teta) - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
     1) * (p1^(teta) * log(p1) + (1 - p2)^(teta) * log((1 - p2)))) - 
     (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * (log((p1^(teta) + 
         (1 - p2)^(teta) - 1)) * (1/teta^2))) * ((1/teta) * (p1^((teta) - 
     1) * (teta))) + (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 
     1) * ((1/teta) * (p1^((teta) - 1) * log(p1) * (teta) + p1^((teta) - 
     1)) - 1/teta^2 * (p1^((teta) - 1) * (teta)))))*(-(exp(teta.st)))
 
 
 c.copula2.be2th <- (((p1^(teta) + (1 - p2)^(teta) - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
     1) * (p1^(teta) * log(p1) + (1 - p2)^(teta) * log((1 - p2)))) - 
     (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * (log((p1^(teta) + 
         (1 - p2)^(teta) - 1)) * (1/teta^2))) * ((1/teta) * ((1 - 
     p2)^((teta) - 1) * (teta))) + (p1^(teta) + (1 - p2)^(teta) - 
     1)^((1/teta) - 1) * ((1/teta) * ((1 - p2)^((teta) - 1) * 
     log((1 - p2)) * (teta) + (1 - p2)^((teta) - 1)) - 1/teta^2 * 
     ((1 - p2)^((teta) - 1) * (teta))))*(-(exp(teta.st)))
 
 bit1.th2 <--(((p1^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
     epsilon)) - 1)^(1/(-(exp(teta.st) + epsilon))) * (log((p1^(-(exp(teta.st) + 
     epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)) * 
     (exp(teta.st)/(-(exp(teta.st) + epsilon))^2)) - (p1^(-(exp(teta.st) + 
     epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)^((1/(-(exp(teta.st) + 
     epsilon))) - 1) * ((1/(-(exp(teta.st) + epsilon))) * ((1 - 
     p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * exp(teta.st)) + 
     p1^(-(exp(teta.st) + epsilon)) * (log(p1) * exp(teta.st))))) * 
     (log((p1^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
         epsilon)) - 1)) * (exp(teta.st)/(-(exp(teta.st) + epsilon))^2)) + 
     (p1^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
         epsilon)) - 1)^(1/(-(exp(teta.st) + epsilon))) * (log((p1^(-(exp(teta.st) + 
         epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)) * 
         (exp(teta.st)/(-(exp(teta.st) + epsilon))^2 + exp(teta.st) * 
             (2 * (exp(teta.st) * (-(exp(teta.st) + epsilon))))/((-(exp(teta.st) + 
             epsilon))^2)^2) - ((1 - p2)^(-(exp(teta.st) + epsilon)) * 
         (log((1 - p2)) * exp(teta.st)) + p1^(-(exp(teta.st) + 
         epsilon)) * (log(p1) * exp(teta.st)))/(p1^(-(exp(teta.st) + 
         epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1) * 
         (exp(teta.st)/(-(exp(teta.st) + epsilon))^2)) - (((p1^(-(exp(teta.st) + 
     epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)^((1/(-(exp(teta.st) + 
     epsilon))) - 1) * (log((p1^(-(exp(teta.st) + epsilon)) + 
     (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(-(exp(teta.st) + 
     epsilon))^2)) - (p1^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
     epsilon)) - 1)^(((1/(-(exp(teta.st) + epsilon))) - 1) - 1) * 
     (((1/(-(exp(teta.st) + epsilon))) - 1) * ((1 - p2)^(-(exp(teta.st) + 
         epsilon)) * (log((1 - p2)) * exp(teta.st)) + p1^(-(exp(teta.st) + 
         epsilon)) * (log(p1) * exp(teta.st))))) * ((1/(-(exp(teta.st) + 
     epsilon))) * ((1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - 
     p2)) * exp(teta.st)) + p1^(-(exp(teta.st) + epsilon)) * (log(p1) * 
     exp(teta.st)))) + (p1^(-(exp(teta.st) + epsilon)) + (1 - 
     p2)^(-(exp(teta.st) + epsilon)) - 1)^((1/(-(exp(teta.st) + 
     epsilon))) - 1) * (exp(teta.st)/(-(exp(teta.st) + epsilon))^2 * 
     ((1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * 
         exp(teta.st)) + p1^(-(exp(teta.st) + epsilon)) * (log(p1) * 
         exp(teta.st))) + (1/(-(exp(teta.st) + epsilon))) * ((1 - 
     p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * exp(teta.st)) - 
     (1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * exp(teta.st)) * 
         (log((1 - p2)) * exp(teta.st)) + (p1^(-(exp(teta.st) + 
     epsilon)) * (log(p1) * exp(teta.st)) - p1^(-(exp(teta.st) + 
     epsilon)) * (log(p1) * exp(teta.st)) * (log(p1) * exp(teta.st)))))))
 
 



}


if(BivD=="F"){


  c.copula.be1 <- -(-1/teta * (1/(1 - exp(-teta)) * (exp(-teta * p1) * teta * (1 - 
    exp(-teta * p2)))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))))


 
  c.copula.be2 <- -(-1/teta * (1/(1 - exp(-teta)) * ((1 - exp(-teta * p1)) * (exp(-teta * 
    p2) * teta))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))))



  c.copula.theta <- 1/teta^2 * log(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2)))) + -1/teta * ((1/(1 - 
    exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - exp(-teta)) * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2)))))


 c.copula2.be1 <-  -1/teta * (1/(1 - exp(-teta)) * (exp(-teta * p1) * teta * teta * 
    (1 - exp(-teta * p2)))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) - 1/(1 - 
    exp(-teta)) * (exp(-teta * p1) * teta * (1 - exp(-teta * 
    p2))) * (1/(1 - exp(-teta)) * (exp(-teta * p1) * teta * (1 - 
    exp(-teta * p2))))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))^2)
                  
 c.copula2.be2 <- -1/teta * (1/(1 - exp(-teta)) * ((1 - exp(-teta * p1)) * (exp(-teta * 
    p2) * teta * teta))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) - 1/(1 - 
    exp(-teta)) * ((1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    teta)) * (1/(1 - exp(-teta)) * ((1 - exp(-teta * p1)) * (exp(-teta * 
    p2) * teta)))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))^2)


c.copula2.be1be2 <- -(-1/teta * (1/(1 - exp(-teta)) * (exp(-teta * p1) * teta * (exp(-teta * 
    p2) * teta))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2)))) + 1/(1 - exp(-teta)) * 
    (exp(-teta * p1) * teta * (1 - exp(-teta * p2))) * (1/(1 - 
    exp(-teta)) * ((1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    teta)))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * 
    p1)) * (1 - exp(-teta * p2))))^2))

c.copula2.be1th <- -(1/teta^2 * (1/(1 - exp(-teta)) * (exp(-teta * p1) * teta * 
    (1 - exp(-teta * p2)))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))) + -1/teta * 
    ((1/(1 - exp(-teta)) * ((exp(-teta * p1) - exp(-teta * p1) * 
        p1 * teta) * (1 - exp(-teta * p2)) + exp(-teta * p1) * 
        teta * (exp(-teta * p2) * p2)) - exp(-teta)/(1 - exp(-teta))^2 * 
        (exp(-teta * p1) * teta * (1 - exp(-teta * p2))))/(1/(1 - 
        exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * 
        (1 - exp(-teta * p2)))) - 1/(1 - exp(-teta)) * (exp(-teta * 
        p1) * teta * (1 - exp(-teta * p2))) * (1/(1 - exp(-teta)) * 
        (exp(-teta) - (exp(-teta * p1) * p1 * (1 - exp(-teta * 
            p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
            p2))) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - exp(-teta)) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - 
        exp(-teta)) * ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * 
        (1 - exp(-teta * p2))))^2))

c.copula2.be2th <- -(1/teta^2 * (1/(1 - exp(-teta)) * ((1 - exp(-teta * p1)) * (exp(-teta * 
    p2) * teta))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))) + -1/teta * ((1/(1 - 
    exp(-teta)) * (exp(-teta * p1) * p1 * (exp(-teta * p2) * 
    teta) + (1 - exp(-teta * p1)) * (exp(-teta * p2) - exp(-teta * 
    p2) * p2 * teta)) - exp(-teta)/(1 - exp(-teta))^2 * ((1 - 
    exp(-teta * p1)) * (exp(-teta * p2) * teta)))/(1/(1 - exp(-teta)) * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2)))) - 1/(1 - exp(-teta)) * ((1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * teta)) * (1/(1 - exp(-teta)) * (exp(-teta) - 
    (exp(-teta * p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * 
        p1)) * (exp(-teta * p2) * p2))) - exp(-teta)/(1 - exp(-teta))^2 * 
    ((1 - exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))/(1/(1 - exp(-teta)) * ((1 - exp(-teta)) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))^2))


bit1.th2 <-1/teta^2 * ((1/(1 - exp(-teta)) * (exp(-teta) - (exp(-teta * 
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


if(BivD=="G0"){

  c.copula.be1 <- exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
    (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p1))^(teta - 
    1) * (teta * (1/p1)))))



 
  c.copula.be2 <- exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
    (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p2))^(teta - 
    1) * (teta * (1/p2)))))





  c.copula.theta <-  (-(exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
    (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p1))^teta * 
    log((-log(p1))) + (-log(p2))^teta * log((-log(p2))))) - ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta) * (log(((-log(p1))^teta + (-log(p2))^teta)) * 
    (1/teta^2)))))*exp(teta.st)
    
c.copula2.be1 <- exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
    (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p1))^(teta - 
    1) * (teta * (1/p1))))) * (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 
    1) * ((1/teta) * ((-log(p1))^(teta - 1) * (teta * (1/p1))))) - 
    exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
        (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p1))^(teta - 
        1) * (teta * (1/p1^2)) + (-log(p1))^((teta - 1) - 1) * 
        ((teta - 1) * (1/p1)) * (teta * (1/p1)))) + ((-log(p1))^teta + 
        (-log(p2))^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 
        1) * ((-log(p1))^(teta - 1) * (teta * (1/p1)))) * ((1/teta) * 
        ((-log(p1))^(teta - 1) * (teta * (1/p1)))))


                  
c.copula2.be2 <- exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
    (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p2))^(teta - 
    1) * (teta * (1/p2))))) * (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 
    1) * ((1/teta) * ((-log(p2))^(teta - 1) * (teta * (1/p2))))) - 
    exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
        (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p2))^(teta - 
        1) * (teta * (1/p2^2)) + (-log(p2))^((teta - 1) - 1) * 
        ((teta - 1) * (1/p2)) * (teta * (1/p2)))) + ((-log(p1))^teta + 
        (-log(p2))^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 
        1) * ((-log(p2))^(teta - 1) * (teta * (1/p2)))) * ((1/teta) * 
        ((-log(p2))^(teta - 1) * (teta * (1/p2)))))


c.copula2.be1be2 <- exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
    (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p2))^(teta - 
    1) * (teta * (1/p2))))) * (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 
    1) * ((1/teta) * ((-log(p1))^(teta - 1) * (teta * (1/p1))))) - 
    exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + 
        (-log(p2))^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 
        1) * ((-log(p2))^(teta - 1) * (teta * (1/p2)))) * ((1/teta) * 
        ((-log(p1))^(teta - 1) * (teta * (1/p1)))))



c.copula2.be1th <-(exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * ((((-log(p1))^teta + 
    (-log(p2))^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 1) * 
    ((-log(p1))^teta * log((-log(p1))) + (-log(p2))^teta * log((-log(p2))))) - 
    ((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 1) * (log(((-log(p1))^teta + 
        (-log(p2))^teta)) * (1/teta^2))) * ((1/teta) * ((-log(p1))^(teta - 
    1) * (teta * (1/p1)))) + ((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 
    1) * ((1/teta) * ((-log(p1))^(teta - 1) * log((-log(p1))) * 
    (teta * (1/p1)) + (-log(p1))^(teta - 1) * (1/p1)) - 1/teta^2 * 
    ((-log(p1))^(teta - 1) * (teta * (1/p1))))) - exp(-((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 
    1) * ((1/teta) * ((-log(p1))^teta * log((-log(p1))) + (-log(p2))^teta * 
    log((-log(p2))))) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
    (log(((-log(p1))^teta + (-log(p2))^teta)) * (1/teta^2))) * 
    (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * 
        ((-log(p1))^(teta - 1) * (teta * (1/p1))))))*exp(teta.st)

c.copula2.be2th <-(exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * ((((-log(p1))^teta + 
    (-log(p2))^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 1) * 
    ((-log(p1))^teta * log((-log(p1))) + (-log(p2))^teta * log((-log(p2))))) - 
    ((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 1) * (log(((-log(p1))^teta + 
        (-log(p2))^teta)) * (1/teta^2))) * ((1/teta) * ((-log(p2))^(teta - 
    1) * (teta * (1/p2)))) + ((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 
    1) * ((1/teta) * ((-log(p2))^(teta - 1) * log((-log(p2))) * 
    (teta * (1/p2)) + (-log(p2))^(teta - 1) * (1/p2)) - 1/teta^2 * 
    ((-log(p2))^(teta - 1) * (teta * (1/p2))))) - exp(-((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 
    1) * ((1/teta) * ((-log(p1))^teta * log((-log(p1))) + (-log(p2))^teta * 
    log((-log(p2))))) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
    (log(((-log(p1))^teta + (-log(p2))^teta)) * (1/teta^2))) * 
    (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * 
        ((-log(p2))^(teta - 1) * (teta * (1/p2))))))*exp(teta.st)

bit1.th2 <- -(exp(-((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))^(1/((exp(teta.st) + 1)))) * ((((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))^(((1/((exp(teta.st) + 
    1))) - 1) - 1) * (((1/((exp(teta.st) + 1))) - 1) * ((-log(p1))^((exp(teta.st) + 
    1)) * (log((-log(p1))) * exp(teta.st)) + (-log(p2))^((exp(teta.st) + 
    1)) * (log((-log(p2))) * exp(teta.st)))) - ((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))^((1/((exp(teta.st) + 
    1))) - 1) * (log(((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))) * (exp(teta.st)/((exp(teta.st) + 1))^2))) * ((1/((exp(teta.st) + 
    1))) * ((-log(p1))^((exp(teta.st) + 1)) * (log((-log(p1))) * 
    exp(teta.st)) + (-log(p2))^((exp(teta.st) + 1)) * (log((-log(p2))) * 
    exp(teta.st)))) + ((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))^((1/((exp(teta.st) + 1))) - 1) * ((1/((exp(teta.st) + 
    1))) * ((-log(p1))^((exp(teta.st) + 1)) * (log((-log(p1))) * 
    exp(teta.st)) * (log((-log(p1))) * exp(teta.st)) + (-log(p1))^((exp(teta.st) + 
    1)) * (log((-log(p1))) * exp(teta.st)) + ((-log(p2))^((exp(teta.st) + 
    1)) * (log((-log(p2))) * exp(teta.st)) * (log((-log(p2))) * 
    exp(teta.st)) + (-log(p2))^((exp(teta.st) + 1)) * (log((-log(p2))) * 
    exp(teta.st)))) - exp(teta.st)/((exp(teta.st) + 1))^2 * ((-log(p1))^((exp(teta.st) + 
    1)) * (log((-log(p1))) * exp(teta.st)) + (-log(p2))^((exp(teta.st) + 
    1)) * (log((-log(p2))) * exp(teta.st)))) - ((((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))^((1/((exp(teta.st) + 
    1))) - 1) * ((1/((exp(teta.st) + 1))) * ((-log(p1))^((exp(teta.st) + 
    1)) * (log((-log(p1))) * exp(teta.st)) + (-log(p2))^((exp(teta.st) + 
    1)) * (log((-log(p2))) * exp(teta.st)))) - ((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))^(1/((exp(teta.st) + 
    1))) * (log(((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))) * (exp(teta.st)/((exp(teta.st) + 1))^2))) * (log(((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))) * (exp(teta.st)/((exp(teta.st) + 
    1))^2)) + ((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))^(1/((exp(teta.st) + 1))) * (((-log(p1))^((exp(teta.st) + 
    1)) * (log((-log(p1))) * exp(teta.st)) + (-log(p2))^((exp(teta.st) + 
    1)) * (log((-log(p2))) * exp(teta.st)))/((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1))) * (exp(teta.st)/((exp(teta.st) + 
    1))^2) + log(((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))) * (exp(teta.st)/((exp(teta.st) + 1))^2 - exp(teta.st) * 
    (2 * (exp(teta.st) * ((exp(teta.st) + 1))))/(((exp(teta.st) + 
    1))^2)^2)))) - exp(-((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))^(1/((exp(teta.st) + 1)))) * (((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))^((1/((exp(teta.st) + 
    1))) - 1) * ((1/((exp(teta.st) + 1))) * ((-log(p1))^((exp(teta.st) + 
    1)) * (log((-log(p1))) * exp(teta.st)) + (-log(p2))^((exp(teta.st) + 
    1)) * (log((-log(p2))) * exp(teta.st)))) - ((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))^(1/((exp(teta.st) + 
    1))) * (log(((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))) * (exp(teta.st)/((exp(teta.st) + 1))^2))) * (((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))^((1/((exp(teta.st) + 
    1))) - 1) * ((1/((exp(teta.st) + 1))) * ((-log(p1))^((exp(teta.st) + 
    1)) * (log((-log(p1))) * exp(teta.st)) + (-log(p2))^((exp(teta.st) + 
    1)) * (log((-log(p2))) * exp(teta.st)))) - ((-log(p1))^((exp(teta.st) + 
    1)) + (-log(p2))^((exp(teta.st) + 1)))^(1/((exp(teta.st) + 
    1))) * (log(((-log(p1))^((exp(teta.st) + 1)) + (-log(p2))^((exp(teta.st) + 
    1)))) * (exp(teta.st)/((exp(teta.st) + 1))^2))))



}



if(BivD=="G90"){


  c.copula.be1 <- exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(1 - p1))^((-teta) - 1) * ((-teta) * 
        (1/(1 - p1))))))

 
  c.copula.be2 <- 1 - exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p2))^((-teta) - 1) * ((-teta) * 
        (1/p2)))))





  c.copula.theta <-    (exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta) * 
        (log(((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))) * 
            (1/teta^2)) - ((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p2))^(-teta) * log((-log(p2))) + 
        (-log(1 - p1))^(-teta) * log((-log(1 - p1)))))))*(-exp(teta.st))
    

c.copula2.be1 <- exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(1 - p1))^((-teta) - 
        1) * ((-teta) * (1/(1 - p1))))) * ((-1/teta) * ((-log(1 - 
        p1))^((-teta) - 1) * ((-teta) * (1/(1 - p1))))) + ((-log(1 - 
        p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 1) * 
        ((-1/teta) * ((-log(1 - p1))^(((-teta) - 1) - 1) * (((-teta) - 
            1) * (1/(1 - p1))) * ((-teta) * (1/(1 - p1))) + (-log(1 - 
            p1))^((-teta) - 1) * ((-teta) * (1/(1 - p1)^2))))) - 
    exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
        (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
            1) * ((-1/teta) * ((-log(1 - p1))^((-teta) - 1) * 
            ((-teta) * (1/(1 - p1)))))) * (((-log(1 - p1))^(-teta) + 
        (-log(p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((-log(1 - 
        p1))^((-teta) - 1) * ((-teta) * (1/(1 - p1))))))




                  
c.copula2.be2 <- -(exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p2))^((-teta) - 1) * ((-teta) * 
        (1/p2))))) * (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((-log(p2))^((-teta) - 1) * ((-teta) * 
    (1/p2))))) - exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p2))^((-teta) - 1) * ((-teta) * 
        (1/p2^2)) + (-log(p2))^(((-teta) - 1) - 1) * (((-teta) - 
        1) * (1/p2)) * ((-teta) * (1/p2)))) + ((-log(1 - p1))^(-teta) + 
        (-log(p2))^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
        1) * ((-log(p2))^((-teta) - 1) * ((-teta) * (1/p2)))) * 
        ((-1/teta) * ((-log(p2))^((-teta) - 1) * ((-teta) * (1/p2))))))
    
    
c.copula2.be1be2 <- exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p2))^((-teta) - 1) * ((-teta) * 
        (1/p2))))) * (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((-log(1 - p1))^((-teta) - 1) * ((-teta) * 
    (1/(1 - p1)))))) - exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(p2))^((-teta) - 1) * 
        ((-teta) * (1/p2)))) * ((-1/teta) * ((-log(1 - p1))^((-teta) - 
        1) * ((-teta) * (1/(1 - p1))))))

c.copula2.be1th <-(exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    ((((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * (log(((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))) * 
        (1/teta^2)) - ((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(p2))^(-teta) * log((-log(p2))) + 
        (-log(1 - p1))^(-teta) * log((-log(1 - p1)))))) * ((-1/teta) * 
        ((-log(1 - p1))^((-teta) - 1) * ((-teta) * (1/(1 - p1))))) + 
        ((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
            1) * (1/teta^2 * ((-log(1 - p1))^((-teta) - 1) * 
            ((-teta) * (1/(1 - p1)))) - (-1/teta) * ((-log(1 - 
            p1))^((-teta) - 1) * (1/(1 - p1)) + (-log(1 - p1))^((-teta) - 
            1) * log((-log(1 - p1))) * ((-teta) * (1/(1 - p1)))))) - 
    exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
        (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta) * 
            (log(((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))) * 
                (1/teta^2)) - ((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
            1) * ((-1/teta) * ((-log(p2))^(-teta) * log((-log(p2))) + 
            (-log(1 - p1))^(-teta) * log((-log(1 - p1)))))) * 
        (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
            1) * ((-1/teta) * ((-log(1 - p1))^((-teta) - 1) * 
            ((-teta) * (1/(1 - p1)))))))*(-exp(teta.st))


c.copula2.be2th <- (-(exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    ((((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * (log(((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))) * 
        (1/teta^2)) - ((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(p2))^(-teta) * log((-log(p2))) + 
        (-log(1 - p1))^(-teta) * log((-log(1 - p1)))))) * ((-1/teta) * 
        ((-log(p2))^((-teta) - 1) * ((-teta) * (1/p2)))) + ((-log(1 - 
        p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 1) * 
        (1/teta^2 * ((-log(p2))^((-teta) - 1) * ((-teta) * (1/p2))) - 
            (-1/teta) * ((-log(p2))^((-teta) - 1) * (1/p2) + 
                (-log(p2))^((-teta) - 1) * log((-log(p2))) * 
                  ((-teta) * (1/p2))))) - exp(-((-log(1 - p1))^(-teta) + 
    (-log(p2))^(-teta))^(-1/teta)) * (((-log(1 - p1))^(-teta) + 
    (-log(p2))^(-teta))^(-1/teta) * (log(((-log(1 - p1))^(-teta) + 
    (-log(p2))^(-teta))) * (1/teta^2)) - ((-log(1 - p1))^(-teta) + 
    (-log(p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((-log(p2))^(-teta) * 
    log((-log(p2))) + (-log(1 - p1))^(-teta) * log((-log(1 - 
    p1)))))) * (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((-log(p2))^((-teta) - 1) * ((-teta) * 
    (1/p2)))))))*(-exp(teta.st))


bit1.th2 <-  exp(-((-log(1 - p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
    1))^(1/(exp(teta.st) + 1))) * ((((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(p2))^(exp(teta.st) + 1))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * ((-log(1 - 
    p1))^(exp(teta.st) + 1) * (log((-log(1 - p1))) * exp(teta.st)) + 
    (-log(p2))^(exp(teta.st) + 1) * (log((-log(p2))) * exp(teta.st)))) - 
    ((-log(1 - p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))^((1/(exp(teta.st) + 1)) - 1) * (log(((-log(1 - p1))^(exp(teta.st) + 
        1) + (-log(p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2))) * ((1/(exp(teta.st) + 1)) * ((-log(1 - p1))^(exp(teta.st) + 
    1) * (log((-log(1 - p1))) * exp(teta.st)) + (-log(p2))^(exp(teta.st) + 
    1) * (log((-log(p2))) * exp(teta.st)))) + ((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * ((-log(1 - p1))^(exp(teta.st) + 
    1) * (log((-log(1 - p1))) * exp(teta.st)) * (log((-log(1 - 
    p1))) * exp(teta.st)) + (-log(1 - p1))^(exp(teta.st) + 1) * 
    (log((-log(1 - p1))) * exp(teta.st)) + ((-log(p2))^(exp(teta.st) + 
    1) * (log((-log(p2))) * exp(teta.st)) * (log((-log(p2))) * 
    exp(teta.st)) + (-log(p2))^(exp(teta.st) + 1) * (log((-log(p2))) * 
    exp(teta.st)))) - exp(teta.st)/(exp(teta.st) + 1)^2 * ((-log(1 - 
    p1))^(exp(teta.st) + 1) * (log((-log(1 - p1))) * exp(teta.st)) + 
    (-log(p2))^(exp(teta.st) + 1) * (log((-log(p2))) * exp(teta.st)))) - 
    ((((-log(1 - p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 
        1)) * ((-log(1 - p1))^(exp(teta.st) + 1) * (log((-log(1 - 
        p1))) * exp(teta.st)) + (-log(p2))^(exp(teta.st) + 1) * 
        (log((-log(p2))) * exp(teta.st)))) - ((-log(1 - p1))^(exp(teta.st) + 
        1) + (-log(p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 
        1)) * (log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))) * (exp(teta.st)/(exp(teta.st) + 1)^2))) * (log(((-log(1 - 
        p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))) * (exp(teta.st)/(exp(teta.st) + 1)^2)) + ((-log(1 - 
        p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))^(1/(exp(teta.st) + 1)) * (((-log(1 - p1))^(exp(teta.st) + 
        1) * (log((-log(1 - p1))) * exp(teta.st)) + (-log(p2))^(exp(teta.st) + 
        1) * (log((-log(p2))) * exp(teta.st)))/((-log(1 - p1))^(exp(teta.st) + 
        1) + (-log(p2))^(exp(teta.st) + 1)) * (exp(teta.st)/(exp(teta.st) + 
        1)^2) + log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))) * (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * 
        (2 * (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
        1)^2)^2)))) - exp(-((-log(1 - p1))^(exp(teta.st) + 1) + 
    (-log(p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 1))) * 
    (((-log(1 - p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 
        1)) * ((-log(1 - p1))^(exp(teta.st) + 1) * (log((-log(1 - 
        p1))) * exp(teta.st)) + (-log(p2))^(exp(teta.st) + 1) * 
        (log((-log(p2))) * exp(teta.st)))) - ((-log(1 - p1))^(exp(teta.st) + 
        1) + (-log(p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 
        1)) * (log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))) * (exp(teta.st)/(exp(teta.st) + 1)^2))) * (((-log(1 - 
    p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * ((-log(1 - p1))^(exp(teta.st) + 
    1) * (log((-log(1 - p1))) * exp(teta.st)) + (-log(p2))^(exp(teta.st) + 
    1) * (log((-log(p2))) * exp(teta.st)))) - ((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 1)) * 
    (log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(p2))^(exp(teta.st) + 
        1))) * (exp(teta.st)/(exp(teta.st) + 1)^2)))



}


if(BivD=="G180"){



 
  c.copula.be1 <- 1 - exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
        1) * ((1/teta) * ((-log(1 - p1))^(teta - 1) * (teta * 
        (1/(1 - p1))))))



 
  c.copula.be2 <- 1 - exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
        1) * ((1/teta) * ((-log(1 - p2))^(teta - 1) * (teta * 
        (1/(1 - p2))))))





  c.copula.theta <-   (-(exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
        1) * ((1/teta) * ((-log(1 - p1))^teta * log((-log(1 - 
        p1))) + (-log(1 - p2))^teta * log((-log(1 - p2))))) - 
        ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta) * 
            (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * 
                (1/teta^2)))))*(exp(teta.st))

  
c.copula2.be1 <- -(exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * ((-log(1 - p1))^(teta - 1) * 
        (teta * (1/(1 - p1))))) * ((1/teta) * ((-log(1 - p1))^(teta - 
        1) * (teta * (1/(1 - p1))))) + ((-log(1 - p1))^teta + 
        (-log(1 - p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - 
        p1))^((teta - 1) - 1) * ((teta - 1) * (1/(1 - p1))) * 
        (teta * (1/(1 - p1))) + (-log(1 - p1))^(teta - 1) * (teta * 
        (1/(1 - p1)^2))))) - exp(-((-log(1 - p1))^teta + (-log(1 - 
    p2))^teta)^(1/teta)) * (((-log(1 - p1))^teta + (-log(1 - 
    p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - p1))^(teta - 
    1) * (teta * (1/(1 - p1)))))) * (((-log(1 - p1))^teta + (-log(1 - 
    p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - p1))^(teta - 
    1) * (teta * (1/(1 - p1)))))))




                  
c.copula2.be2 <- -(exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * ((-log(1 - p2))^(teta - 1) * 
        (teta * (1/(1 - p2))))) * ((1/teta) * ((-log(1 - p2))^(teta - 
        1) * (teta * (1/(1 - p2))))) + ((-log(1 - p1))^teta + 
        (-log(1 - p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - 
        p2))^((teta - 1) - 1) * ((teta - 1) * (1/(1 - p2))) * 
        (teta * (1/(1 - p2))) + (-log(1 - p2))^(teta - 1) * (teta * 
        (1/(1 - p2)^2))))) - exp(-((-log(1 - p1))^teta + (-log(1 - 
    p2))^teta)^(1/teta)) * (((-log(1 - p1))^teta + (-log(1 - 
    p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - p2))^(teta - 
    1) * (teta * (1/(1 - p2)))))) * (((-log(1 - p1))^teta + (-log(1 - 
    p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - p2))^(teta - 
    1) * (teta * (1/(1 - p2)))))))

c.copula2.be1be2 <- -(exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * ((-log(1 - p2))^(teta - 1) * 
        (teta * (1/(1 - p2))))) * ((1/teta) * ((-log(1 - p1))^(teta - 
        1) * (teta * (1/(1 - p1)))))) - exp(-((-log(1 - p1))^teta + 
    (-log(1 - p2))^teta)^(1/teta)) * (((-log(1 - p1))^teta + 
    (-log(1 - p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - 
    p2))^(teta - 1) * (teta * (1/(1 - p2)))))) * (((-log(1 - 
    p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 1) * ((1/teta) * 
    ((-log(1 - p1))^(teta - 1) * (teta * (1/(1 - p1)))))))


c.copula2.be1th <-(-(exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    ((((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * ((-log(1 - p1))^teta * log((-log(1 - 
        p1))) + (-log(1 - p2))^teta * log((-log(1 - p2))))) - 
        ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
            1) * (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * 
            (1/teta^2))) * ((1/teta) * ((-log(1 - p1))^(teta - 
        1) * (teta * (1/(1 - p1))))) + ((-log(1 - p1))^teta + 
        (-log(1 - p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - 
        p1))^(teta - 1) * log((-log(1 - p1))) * (teta * (1/(1 - 
        p1))) + (-log(1 - p1))^(teta - 1) * (1/(1 - p1))) - 1/teta^2 * 
        ((-log(1 - p1))^(teta - 1) * (teta * (1/(1 - p1)))))) - 
    exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
        (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
            1) * ((1/teta) * ((-log(1 - p1))^teta * log((-log(1 - 
            p1))) + (-log(1 - p2))^teta * log((-log(1 - p2))))) - 
            ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta) * 
                (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * 
                  (1/teta^2))) * (((-log(1 - p1))^teta + (-log(1 - 
        p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - p1))^(teta - 
        1) * (teta * (1/(1 - p1))))))))*exp(teta.st)


c.copula2.be2th <- (-(exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    ((((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * ((-log(1 - p1))^teta * log((-log(1 - 
        p1))) + (-log(1 - p2))^teta * log((-log(1 - p2))))) - 
        ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
            1) * (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * 
            (1/teta^2))) * ((1/teta) * ((-log(1 - p2))^(teta - 
        1) * (teta * (1/(1 - p2))))) + ((-log(1 - p1))^teta + 
        (-log(1 - p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - 
        p2))^(teta - 1) * log((-log(1 - p2))) * (teta * (1/(1 - 
        p2))) + (-log(1 - p2))^(teta - 1) * (1/(1 - p2))) - 1/teta^2 * 
        ((-log(1 - p2))^(teta - 1) * (teta * (1/(1 - p2)))))) - 
    exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
        (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
            1) * ((1/teta) * ((-log(1 - p1))^teta * log((-log(1 - 
            p1))) + (-log(1 - p2))^teta * log((-log(1 - p2))))) - 
            ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta) * 
                (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * 
                  (1/teta^2))) * (((-log(1 - p1))^teta + (-log(1 - 
        p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - p2))^(teta - 
        1) * (teta * (1/(1 - p2))))))))*exp(teta.st)


bit1.th2 <- -(exp(-((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))^(1/(exp(teta.st) + 1))) * ((((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * ((-log(1 - 
    p1))^(exp(teta.st) + 1) * (log((-log(1 - p1))) * exp(teta.st)) + 
    (-log(1 - p2))^(exp(teta.st) + 1) * (log((-log(1 - p2))) * 
        exp(teta.st)))) - ((-log(1 - p1))^(exp(teta.st) + 1) + 
    (-log(1 - p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 1)) - 
    1) * (log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - 
    p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2))) * ((1/(exp(teta.st) + 1)) * ((-log(1 - p1))^(exp(teta.st) + 
    1) * (log((-log(1 - p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))) + ((-log(1 - 
    p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
    ((-log(1 - p1))^(exp(teta.st) + 1) * (log((-log(1 - p1))) * 
        exp(teta.st)) * (log((-log(1 - p1))) * exp(teta.st)) + 
        (-log(1 - p1))^(exp(teta.st) + 1) * (log((-log(1 - p1))) * 
            exp(teta.st)) + ((-log(1 - p2))^(exp(teta.st) + 1) * 
        (log((-log(1 - p2))) * exp(teta.st)) * (log((-log(1 - 
        p2))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
        1) * (log((-log(1 - p2))) * exp(teta.st)))) - exp(teta.st)/(exp(teta.st) + 
    1)^2 * ((-log(1 - p1))^(exp(teta.st) + 1) * (log((-log(1 - 
    p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 1) * 
    (log((-log(1 - p2))) * exp(teta.st)))) - ((((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * ((-log(1 - p1))^(exp(teta.st) + 
    1) * (log((-log(1 - p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))) - ((-log(1 - 
    p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))^(1/(exp(teta.st) + 1)) * (log(((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2))) * (log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - 
    p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + ((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))^(1/(exp(teta.st) + 1)) * (((-log(1 - p1))^(exp(teta.st) + 
    1) * (log((-log(1 - p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))/((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1)) * (exp(teta.st)/(exp(teta.st) + 
    1)^2) + log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - 
    p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
    1)))/((exp(teta.st) + 1)^2)^2)))) - exp(-((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 
    1))) * (((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
    ((-log(1 - p1))^(exp(teta.st) + 1) * (log((-log(1 - p1))) * 
        exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 1) * (log((-log(1 - 
        p2))) * exp(teta.st)))) - ((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 
    1)) * (log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - 
    p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2))) * (((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - 
    p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 
    1)) * ((-log(1 - p1))^(exp(teta.st) + 1) * (log((-log(1 - 
    p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 1) * 
    (log((-log(1 - p2))) * exp(teta.st)))) - ((-log(1 - p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 
    1)) * (log(((-log(1 - p1))^(exp(teta.st) + 1) + (-log(1 - 
    p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2))))






}


if(BivD=="G270"){


   c.copula.be1 <- 1 - exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
     (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
         1) * ((-1/teta) * ((-log(p1))^((-teta) - 1) * ((-teta) * 
         (1/p1)))))
 
 
  
   c.copula.be2 <- exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
     (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
         1) * ((-1/teta) * ((-log(1 - p2))^((-teta) - 1) * ((-teta) * 
         (1/(1 - p2))))))
 
 
 
 
   c.copula.theta <-  (exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
     (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta) * 
         (log(((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))) * 
             (1/teta^2)) - ((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
         1) * ((-1/teta) * ((-log(1 - p2))^(-teta) * log((-log(1 - 
         p2))) + (-log(p1))^(-teta) * log((-log(p1)))))))*(-exp(teta.st))
 
c.copula2.be1 <- -(exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
    (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p1))^((-teta) - 1) * ((-teta) * 
        (1/p1))))) * (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((-log(p1))^((-teta) - 1) * ((-teta) * 
    (1/p1))))) - exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
    (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p1))^((-teta) - 1) * ((-teta) * 
        (1/p1^2)) + (-log(p1))^(((-teta) - 1) - 1) * (((-teta) - 
        1) * (1/p1)) * ((-teta) * (1/p1)))) + ((-log(p1))^(-teta) + 
        (-log(1 - p2))^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
        1) * ((-log(p1))^((-teta) - 1) * ((-teta) * (1/p1)))) * 
        ((-1/teta) * ((-log(p1))^((-teta) - 1) * ((-teta) * (1/p1))))))


                  
c.copula2.be2 <- exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
    (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(1 - p2))^((-teta) - 
        1) * ((-teta) * (1/(1 - p2))))) * ((-1/teta) * ((-log(1 - 
        p2))^((-teta) - 1) * ((-teta) * (1/(1 - p2))))) + ((-log(p1))^(-teta) + 
        (-log(1 - p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * 
        ((-log(1 - p2))^(((-teta) - 1) - 1) * (((-teta) - 1) * 
            (1/(1 - p2))) * ((-teta) * (1/(1 - p2))) + (-log(1 - 
            p2))^((-teta) - 1) * ((-teta) * (1/(1 - p2)^2))))) - 
    exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
        (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
            1) * ((-1/teta) * ((-log(1 - p2))^((-teta) - 1) * 
            ((-teta) * (1/(1 - p2)))))) * (((-log(p1))^(-teta) + 
        (-log(1 - p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * 
        ((-log(1 - p2))^((-teta) - 1) * ((-teta) * (1/(1 - p2))))))



c.copula2.be1be2 <- -(exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
    (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(1 - p2))^((-teta) - 
        1) * ((-teta) * (1/(1 - p2))))) * ((-1/teta) * ((-log(p1))^((-teta) - 
        1) * ((-teta) * (1/p1))))) - exp(-((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^(-1/teta)) * (((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((-log(1 - 
    p2))^((-teta) - 1) * ((-teta) * (1/(1 - p2)))))) * (((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((-log(p1))^((-teta) - 
    1) * ((-teta) * (1/p1))))))

c.copula2.be1th <-(-(exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
    ((((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
        1) * (log(((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))) * 
        (1/teta^2)) - ((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(1 - p2))^(-teta) * 
        log((-log(1 - p2))) + (-log(p1))^(-teta) * log((-log(p1)))))) * 
        ((-1/teta) * ((-log(p1))^((-teta) - 1) * ((-teta) * (1/p1)))) + 
        ((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
            1) * (1/teta^2 * ((-log(p1))^((-teta) - 1) * ((-teta) * 
            (1/p1))) - (-1/teta) * ((-log(p1))^((-teta) - 1) * 
            (1/p1) + (-log(p1))^((-teta) - 1) * log((-log(p1))) * 
            ((-teta) * (1/p1))))) - exp(-((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^(-1/teta)) * (((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^(-1/teta) * (log(((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))) * (1/teta^2)) - ((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((-log(1 - 
    p2))^(-teta) * log((-log(1 - p2))) + (-log(p1))^(-teta) * 
    log((-log(p1)))))) * (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((-log(p1))^((-teta) - 1) * ((-teta) * 
    (1/p1)))))))*(-exp(teta.st))



c.copula2.be2th <- (exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
    ((((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
        1) * (log(((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))) * 
        (1/teta^2)) - ((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(1 - p2))^(-teta) * 
        log((-log(1 - p2))) + (-log(p1))^(-teta) * log((-log(p1)))))) * 
        ((-1/teta) * ((-log(1 - p2))^((-teta) - 1) * ((-teta) * 
            (1/(1 - p2))))) + ((-log(p1))^(-teta) + (-log(1 - 
        p2))^(-teta))^((-1/teta) - 1) * (1/teta^2 * ((-log(1 - 
        p2))^((-teta) - 1) * ((-teta) * (1/(1 - p2)))) - (-1/teta) * 
        ((-log(1 - p2))^((-teta) - 1) * (1/(1 - p2)) + (-log(1 - 
            p2))^((-teta) - 1) * log((-log(1 - p2))) * ((-teta) * 
            (1/(1 - p2)))))) - exp(-((-log(p1))^(-teta) + (-log(1 - 
    p2))^(-teta))^(-1/teta)) * (((-log(p1))^(-teta) + (-log(1 - 
    p2))^(-teta))^(-1/teta) * (log(((-log(p1))^(-teta) + (-log(1 - 
    p2))^(-teta))) * (1/teta^2)) - ((-log(p1))^(-teta) + (-log(1 - 
    p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((-log(1 - p2))^(-teta) * 
    log((-log(1 - p2))) + (-log(p1))^(-teta) * log((-log(p1)))))) * 
    (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(1 - p2))^((-teta) - 1) * ((-teta) * 
        (1/(1 - p2)))))))*(-exp(teta.st))



bit1.th2 <-exp(-((-log(p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))^(1/(exp(teta.st) + 1))) * ((((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * ((-log(p1))^(exp(teta.st) + 
    1) * (log((-log(p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))) - ((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 
    1)) - 1) * (log(((-log(p1))^(exp(teta.st) + 1) + (-log(1 - 
    p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2))) * ((1/(exp(teta.st) + 1)) * ((-log(p1))^(exp(teta.st) + 
    1) * (log((-log(p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))) + ((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * ((-log(p1))^(exp(teta.st) + 
    1) * (log((-log(p1))) * exp(teta.st)) * (log((-log(p1))) * 
    exp(teta.st)) + (-log(p1))^(exp(teta.st) + 1) * (log((-log(p1))) * 
    exp(teta.st)) + ((-log(1 - p2))^(exp(teta.st) + 1) * (log((-log(1 - 
    p2))) * exp(teta.st)) * (log((-log(1 - p2))) * exp(teta.st)) + 
    (-log(1 - p2))^(exp(teta.st) + 1) * (log((-log(1 - p2))) * 
        exp(teta.st)))) - exp(teta.st)/(exp(teta.st) + 1)^2 * 
    ((-log(p1))^(exp(teta.st) + 1) * (log((-log(p1))) * exp(teta.st)) + 
        (-log(1 - p2))^(exp(teta.st) + 1) * (log((-log(1 - p2))) * 
            exp(teta.st)))) - ((((-log(p1))^(exp(teta.st) + 1) + 
    (-log(1 - p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * ((-log(p1))^(exp(teta.st) + 
    1) * (log((-log(p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))) - ((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 
    1)) * (log(((-log(p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))) * (exp(teta.st)/(exp(teta.st) + 1)^2))) * (log(((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + ((-log(p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))^(1/(exp(teta.st) + 1)) * (((-log(p1))^(exp(teta.st) + 
    1) * (log((-log(p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))/((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1)) * (exp(teta.st)/(exp(teta.st) + 
    1)^2) + log(((-log(p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))) * (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * 
    (2 * (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
    1)^2)^2)))) - exp(-((-log(p1))^(exp(teta.st) + 1) + (-log(1 - 
    p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 1))) * (((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * ((-log(p1))^(exp(teta.st) + 
    1) * (log((-log(p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))) - ((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 
    1)) * (log(((-log(p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))) * (exp(teta.st)/(exp(teta.st) + 1)^2))) * (((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * ((-log(p1))^(exp(teta.st) + 
    1) * (log((-log(p1))) * exp(teta.st)) + (-log(1 - p2))^(exp(teta.st) + 
    1) * (log((-log(1 - p2))) * exp(teta.st)))) - ((-log(p1))^(exp(teta.st) + 
    1) + (-log(1 - p2))^(exp(teta.st) + 1))^(1/(exp(teta.st) + 
    1)) * (log(((-log(p1))^(exp(teta.st) + 1) + (-log(1 - p2))^(exp(teta.st) + 
    1))) * (exp(teta.st)/(exp(teta.st) + 1)^2)))




}



if(BivD=="J0"){


  c.copula.be1 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p1)^(teta - 1) * teta - (1 - p1)^(teta - 
    1) * teta * (1 - p2)^teta))


 
  c.copula.be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
    ((1 - p2)^(teta - 1) * teta)))




  c.copula.theta <-    (-(((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p1)^teta * log((1 - p1)) + (1 - p2)^teta * 
    log((1 - p2)) - ((1 - p1)^teta * log((1 - p1)) * (1 - p2)^teta + 
    (1 - p1)^teta * ((1 - p2)^teta * log((1 - p2)))))) - ((1 - 
    p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta) * 
    (log(((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - 
        p2)^teta)) * (1/teta^2))))*exp(teta.st)
  



c.copula2.be1 <- -(((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p1)^((teta - 1) - 1) * (teta - 1) * 
    teta - (1 - p1)^((teta - 1) - 1) * (teta - 1) * teta * (1 - 
    p2)^teta)) + ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 1) * ((1 - 
    p1)^(teta - 1) * teta - (1 - p1)^(teta - 1) * teta * (1 - 
    p2)^teta)) * ((1/teta) * ((1 - p1)^(teta - 1) * teta - (1 - 
    p1)^(teta - 1) * teta * (1 - p2)^teta)))


                  
c.copula2.be2 <- -(((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p2)^((teta - 1) - 1) * (teta - 1) * 
    teta - (1 - p1)^teta * ((1 - p2)^((teta - 1) - 1) * (teta - 
    1) * teta))) + ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 1) * ((1 - 
    p2)^(teta - 1) * teta - (1 - p1)^teta * ((1 - p2)^(teta - 
    1) * teta))) * ((1/teta) * ((1 - p2)^(teta - 1) * teta - 
    (1 - p1)^teta * ((1 - p2)^(teta - 1) * teta))))

c.copula2.be1be2 <- ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 
    1) * ((1/teta) * ((1 - p1)^(teta - 1) * teta * ((1 - p2)^(teta - 
    1) * teta))) - ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 1) * ((1 - 
    p2)^(teta - 1) * teta - (1 - p1)^teta * ((1 - p2)^(teta - 
    1) * teta))) * ((1/teta) * ((1 - p1)^(teta - 1) * teta - 
    (1 - p1)^(teta - 1) * teta * (1 - p2)^teta))


c.copula2.be1th <-((((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - p1)^teta * log((1 - p1)) + 
    (1 - p2)^teta * log((1 - p2)) - ((1 - p1)^teta * log((1 - 
    p1)) * (1 - p2)^teta + (1 - p1)^teta * ((1 - p2)^teta * log((1 - 
    p2)))))) - ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^((1/teta) - 1) * (log(((1 - p1)^teta + (1 - 
    p2)^teta - (1 - p1)^teta * (1 - p2)^teta)) * (1/teta^2))) * 
    ((1/teta) * ((1 - p1)^(teta - 1) * teta - (1 - p1)^(teta - 
        1) * teta * (1 - p2)^teta)) + ((1 - p1)^teta + (1 - p2)^teta - 
    (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 1) * ((1/teta) * 
    ((1 - p1)^(teta - 1) * log((1 - p1)) * teta + (1 - p1)^(teta - 
        1) - (((1 - p1)^(teta - 1) * log((1 - p1)) * teta + (1 - 
        p1)^(teta - 1)) * (1 - p2)^teta + (1 - p1)^(teta - 1) * 
        teta * ((1 - p2)^teta * log((1 - p2))))) - 1/teta^2 * 
    ((1 - p1)^(teta - 1) * teta - (1 - p1)^(teta - 1) * teta * 
        (1 - p2)^teta)))*exp(teta.st)


c.copula2.be2th <- ((((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - p1)^teta * log((1 - p1)) + 
    (1 - p2)^teta * log((1 - p2)) - ((1 - p1)^teta * log((1 - 
    p1)) * (1 - p2)^teta + (1 - p1)^teta * ((1 - p2)^teta * log((1 - 
    p2)))))) - ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^((1/teta) - 1) * (log(((1 - p1)^teta + (1 - 
    p2)^teta - (1 - p1)^teta * (1 - p2)^teta)) * (1/teta^2))) * 
    ((1/teta) * ((1 - p2)^(teta - 1) * teta - (1 - p1)^teta * 
        ((1 - p2)^(teta - 1) * teta))) + ((1 - p1)^teta + (1 - 
    p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^((1/teta) - 1) * 
    ((1/teta) * ((1 - p2)^(teta - 1) * log((1 - p2)) * teta + 
        (1 - p2)^(teta - 1) - ((1 - p1)^teta * log((1 - p1)) * 
        ((1 - p2)^(teta - 1) * teta) + (1 - p1)^teta * ((1 - 
        p2)^(teta - 1) * log((1 - p2)) * teta + (1 - p2)^(teta - 
        1)))) - 1/teta^2 * ((1 - p2)^(teta - 1) * teta - (1 - 
        p1)^teta * ((1 - p2)^(teta - 1) * teta))))*exp(teta.st)


bit1.th2 <- -((((1 - p1)^(exp(teta.st) + 1 + epsilon) + (1 - p2)^(exp(teta.st) + 
    1 + epsilon) - (1 - p1)^(exp(teta.st) + 1 + epsilon) * (1 - 
    p2)^(exp(teta.st) + 1 + epsilon))^(((1/(exp(teta.st) + 1 + 
    epsilon)) - 1) - 1) * (((1/(exp(teta.st) + 1 + epsilon)) - 
    1) * ((1 - p1)^(exp(teta.st) + 1 + epsilon) * (log((1 - p1)) * 
    exp(teta.st)) + (1 - p2)^(exp(teta.st) + 1 + epsilon) * (log((1 - 
    p2)) * exp(teta.st)) - ((1 - p1)^(exp(teta.st) + 1 + epsilon) * 
    (log((1 - p1)) * exp(teta.st)) * (1 - p2)^(exp(teta.st) + 
    1 + epsilon) + (1 - p1)^(exp(teta.st) + 1 + epsilon) * ((1 - 
    p2)^(exp(teta.st) + 1 + epsilon) * (log((1 - p2)) * exp(teta.st)))))) - 
    ((1 - p1)^(exp(teta.st) + 1 + epsilon) + (1 - p2)^(exp(teta.st) + 
        1 + epsilon) - (1 - p1)^(exp(teta.st) + 1 + epsilon) * 
        (1 - p2)^(exp(teta.st) + 1 + epsilon))^((1/(exp(teta.st) + 
        1 + epsilon)) - 1) * (log(((1 - p1)^(exp(teta.st) + 1 + 
        epsilon) + (1 - p2)^(exp(teta.st) + 1 + epsilon) - (1 - 
        p1)^(exp(teta.st) + 1 + epsilon) * (1 - p2)^(exp(teta.st) + 
        1 + epsilon))) * (exp(teta.st)/(exp(teta.st) + 1 + epsilon)^2))) * 
    ((1/(exp(teta.st) + 1 + epsilon)) * ((1 - p1)^(exp(teta.st) + 
        1 + epsilon) * (log((1 - p1)) * exp(teta.st)) + (1 - 
        p2)^(exp(teta.st) + 1 + epsilon) * (log((1 - p2)) * exp(teta.st)) - 
        ((1 - p1)^(exp(teta.st) + 1 + epsilon) * (log((1 - p1)) * 
            exp(teta.st)) * (1 - p2)^(exp(teta.st) + 1 + epsilon) + 
            (1 - p1)^(exp(teta.st) + 1 + epsilon) * ((1 - p2)^(exp(teta.st) + 
                1 + epsilon) * (log((1 - p2)) * exp(teta.st)))))) + 
    ((1 - p1)^(exp(teta.st) + 1 + epsilon) + (1 - p2)^(exp(teta.st) + 
        1 + epsilon) - (1 - p1)^(exp(teta.st) + 1 + epsilon) * 
        (1 - p2)^(exp(teta.st) + 1 + epsilon))^((1/(exp(teta.st) + 
        1 + epsilon)) - 1) * ((1/(exp(teta.st) + 1 + epsilon)) * 
        ((1 - p1)^(exp(teta.st) + 1 + epsilon) * (log((1 - p1)) * 
            exp(teta.st)) * (log((1 - p1)) * exp(teta.st)) + 
            (1 - p1)^(exp(teta.st) + 1 + epsilon) * (log((1 - 
                p1)) * exp(teta.st)) + ((1 - p2)^(exp(teta.st) + 
            1 + epsilon) * (log((1 - p2)) * exp(teta.st)) * (log((1 - 
            p2)) * exp(teta.st)) + (1 - p2)^(exp(teta.st) + 1 + 
            epsilon) * (log((1 - p2)) * exp(teta.st))) - (((1 - 
            p1)^(exp(teta.st) + 1 + epsilon) * (log((1 - p1)) * 
            exp(teta.st)) * (log((1 - p1)) * exp(teta.st)) + 
            (1 - p1)^(exp(teta.st) + 1 + epsilon) * (log((1 - 
                p1)) * exp(teta.st))) * (1 - p2)^(exp(teta.st) + 
            1 + epsilon) + (1 - p1)^(exp(teta.st) + 1 + epsilon) * 
            (log((1 - p1)) * exp(teta.st)) * ((1 - p2)^(exp(teta.st) + 
            1 + epsilon) * (log((1 - p2)) * exp(teta.st))) + 
            ((1 - p1)^(exp(teta.st) + 1 + epsilon) * (log((1 - 
                p1)) * exp(teta.st)) * ((1 - p2)^(exp(teta.st) + 
                1 + epsilon) * (log((1 - p2)) * exp(teta.st))) + 
                (1 - p1)^(exp(teta.st) + 1 + epsilon) * ((1 - 
                  p2)^(exp(teta.st) + 1 + epsilon) * (log((1 - 
                  p2)) * exp(teta.st)) * (log((1 - p2)) * exp(teta.st)) + 
                  (1 - p2)^(exp(teta.st) + 1 + epsilon) * (log((1 - 
                    p2)) * exp(teta.st)))))) - exp(teta.st)/(exp(teta.st) + 
        1 + epsilon)^2 * ((1 - p1)^(exp(teta.st) + 1 + epsilon) * 
        (log((1 - p1)) * exp(teta.st)) + (1 - p2)^(exp(teta.st) + 
        1 + epsilon) * (log((1 - p2)) * exp(teta.st)) - ((1 - 
        p1)^(exp(teta.st) + 1 + epsilon) * (log((1 - p1)) * exp(teta.st)) * 
        (1 - p2)^(exp(teta.st) + 1 + epsilon) + (1 - p1)^(exp(teta.st) + 
        1 + epsilon) * ((1 - p2)^(exp(teta.st) + 1 + epsilon) * 
        (log((1 - p2)) * exp(teta.st)))))) - ((((1 - p1)^(exp(teta.st) + 
    1 + epsilon) + (1 - p2)^(exp(teta.st) + 1 + epsilon) - (1 - 
    p1)^(exp(teta.st) + 1 + epsilon) * (1 - p2)^(exp(teta.st) + 
    1 + epsilon))^((1/(exp(teta.st) + 1 + epsilon)) - 1) * ((1/(exp(teta.st) + 
    1 + epsilon)) * ((1 - p1)^(exp(teta.st) + 1 + epsilon) * 
    (log((1 - p1)) * exp(teta.st)) + (1 - p2)^(exp(teta.st) + 
    1 + epsilon) * (log((1 - p2)) * exp(teta.st)) - ((1 - p1)^(exp(teta.st) + 
    1 + epsilon) * (log((1 - p1)) * exp(teta.st)) * (1 - p2)^(exp(teta.st) + 
    1 + epsilon) + (1 - p1)^(exp(teta.st) + 1 + epsilon) * ((1 - 
    p2)^(exp(teta.st) + 1 + epsilon) * (log((1 - p2)) * exp(teta.st)))))) - 
    ((1 - p1)^(exp(teta.st) + 1 + epsilon) + (1 - p2)^(exp(teta.st) + 
        1 + epsilon) - (1 - p1)^(exp(teta.st) + 1 + epsilon) * 
        (1 - p2)^(exp(teta.st) + 1 + epsilon))^(1/(exp(teta.st) + 
        1 + epsilon)) * (log(((1 - p1)^(exp(teta.st) + 1 + epsilon) + 
        (1 - p2)^(exp(teta.st) + 1 + epsilon) - (1 - p1)^(exp(teta.st) + 
        1 + epsilon) * (1 - p2)^(exp(teta.st) + 1 + epsilon))) * 
        (exp(teta.st)/(exp(teta.st) + 1 + epsilon)^2))) * (log(((1 - 
    p1)^(exp(teta.st) + 1 + epsilon) + (1 - p2)^(exp(teta.st) + 
    1 + epsilon) - (1 - p1)^(exp(teta.st) + 1 + epsilon) * (1 - 
    p2)^(exp(teta.st) + 1 + epsilon))) * (exp(teta.st)/(exp(teta.st) + 
    1 + epsilon)^2)) + ((1 - p1)^(exp(teta.st) + 1 + epsilon) + 
    (1 - p2)^(exp(teta.st) + 1 + epsilon) - (1 - p1)^(exp(teta.st) + 
    1 + epsilon) * (1 - p2)^(exp(teta.st) + 1 + epsilon))^(1/(exp(teta.st) + 
    1 + epsilon)) * (((1 - p1)^(exp(teta.st) + 1 + epsilon) * 
    (log((1 - p1)) * exp(teta.st)) + (1 - p2)^(exp(teta.st) + 
    1 + epsilon) * (log((1 - p2)) * exp(teta.st)) - ((1 - p1)^(exp(teta.st) + 
    1 + epsilon) * (log((1 - p1)) * exp(teta.st)) * (1 - p2)^(exp(teta.st) + 
    1 + epsilon) + (1 - p1)^(exp(teta.st) + 1 + epsilon) * ((1 - 
    p2)^(exp(teta.st) + 1 + epsilon) * (log((1 - p2)) * exp(teta.st)))))/((1 - 
    p1)^(exp(teta.st) + 1 + epsilon) + (1 - p2)^(exp(teta.st) + 
    1 + epsilon) - (1 - p1)^(exp(teta.st) + 1 + epsilon) * (1 - 
    p2)^(exp(teta.st) + 1 + epsilon)) * (exp(teta.st)/(exp(teta.st) + 
    1 + epsilon)^2) + log(((1 - p1)^(exp(teta.st) + 1 + epsilon) + 
    (1 - p2)^(exp(teta.st) + 1 + epsilon) - (1 - p1)^(exp(teta.st) + 
    1 + epsilon) * (1 - p2)^(exp(teta.st) + 1 + epsilon))) * 
    (exp(teta.st)/(exp(teta.st) + 1 + epsilon)^2 - exp(teta.st) * 
        (2 * (exp(teta.st) * (exp(teta.st) + 1 + epsilon)))/((exp(teta.st) + 
        1 + epsilon)^2)^2))))



}

















if(BivD=="J90"){

c.copula.be1 <- (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * (p1^((-teta) - 1) * (-teta) - p1^((-teta) - 
    1) * (-teta) * (1 - p2)^(-teta)))

 
  c.copula.be2 <- 1 - (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - p2)^((-teta) - 1) * (-teta) - p1^(-teta) * 
    ((1 - p2)^((-teta) - 1) * (-teta))))




  c.copula.theta <-    ((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^(-1/teta) * 
    (log((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))) * 
        (1/teta^2)) - (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * 
    (1 - p2)^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((1 - p2)^(-teta) * 
    log((1 - p2)) + p1^(-teta) * log(p1) - (p1^(-teta) * ((1 - 
    p2)^(-teta) * log((1 - p2))) + p1^(-teta) * log(p1) * (1 - 
    p2)^(-teta)))))*(-exp(teta.st))
  
c.copula2.be1 <- (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (p1^((-teta) - 1) * (-teta) - 
    p1^((-teta) - 1) * (-teta) * (1 - p2)^(-teta))) * ((-1/teta) * 
    (p1^((-teta) - 1) * (-teta) - p1^((-teta) - 1) * (-teta) * 
        (1 - p2)^(-teta))) + (p1^(-teta) + (1 - p2)^(-teta) - 
    p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 1) * ((-1/teta) * 
    (p1^(((-teta) - 1) - 1) * ((-teta) - 1) * (-teta) - p1^(((-teta) - 
        1) - 1) * ((-teta) - 1) * (-teta) * (1 - p2)^(-teta)))




                  
c.copula2.be2 <- (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - p2)^(((-teta) - 1) - 1) * ((-teta) - 
    1) * (-teta) - p1^(-teta) * ((1 - p2)^(((-teta) - 1) - 1) * 
    ((-teta) - 1) * (-teta)))) + (p1^(-teta) + (1 - p2)^(-teta) - 
    p1^(-teta) * (1 - p2)^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * ((1 - p2)^((-teta) - 1) * (-teta) - p1^(-teta) * ((1 - 
    p2)^((-teta) - 1) * (-teta)))) * ((-1/teta) * ((1 - p2)^((-teta) - 
    1) * (-teta) - p1^(-teta) * ((1 - p2)^((-teta) - 1) * (-teta))))
    
c.copula2.be1be2 <- (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * (p1^((-teta) - 1) * (-teta) * ((1 - p2)^((-teta) - 
    1) * (-teta)))) - (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * 
    (1 - p2)^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 1) * 
    ((1 - p2)^((-teta) - 1) * (-teta) - p1^(-teta) * ((1 - p2)^((-teta) - 
        1) * (-teta)))) * ((-1/teta) * (p1^((-teta) - 1) * (-teta) - 
    p1^((-teta) - 1) * (-teta) * (1 - p2)^(-teta)))

c.copula2.be1th <-(((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
    1) * (log((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - 
    p2)^(-teta))) * (1/teta^2)) - (p1^(-teta) + (1 - p2)^(-teta) - 
    p1^(-teta) * (1 - p2)^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * ((1 - p2)^(-teta) * log((1 - p2)) + p1^(-teta) * log(p1) - 
    (p1^(-teta) * ((1 - p2)^(-teta) * log((1 - p2))) + p1^(-teta) * 
        log(p1) * (1 - p2)^(-teta))))) * ((-1/teta) * (p1^((-teta) - 
    1) * (-teta) - p1^((-teta) - 1) * (-teta) * (1 - p2)^(-teta))) + 
    (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
        1) * (1/teta^2 * (p1^((-teta) - 1) * (-teta) - p1^((-teta) - 
        1) * (-teta) * (1 - p2)^(-teta)) - (-1/teta) * (p1^((-teta) - 
        1) + p1^((-teta) - 1) * log(p1) * (-teta) - (p1^((-teta) - 
        1) * (-teta) * ((1 - p2)^(-teta) * log((1 - p2))) + (p1^((-teta) - 
        1) + p1^((-teta) - 1) * log(p1) * (-teta)) * (1 - p2)^(-teta)))))*(-exp(teta.st))

c.copula2.be2th <- (-(((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
    1) * (log((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - 
    p2)^(-teta))) * (1/teta^2)) - (p1^(-teta) + (1 - p2)^(-teta) - 
    p1^(-teta) * (1 - p2)^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * ((1 - p2)^(-teta) * log((1 - p2)) + p1^(-teta) * log(p1) - 
    (p1^(-teta) * ((1 - p2)^(-teta) * log((1 - p2))) + p1^(-teta) * 
        log(p1) * (1 - p2)^(-teta))))) * ((-1/teta) * ((1 - p2)^((-teta) - 
    1) * (-teta) - p1^(-teta) * ((1 - p2)^((-teta) - 1) * (-teta)))) + 
    (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
        1) * (1/teta^2 * ((1 - p2)^((-teta) - 1) * (-teta) - 
        p1^(-teta) * ((1 - p2)^((-teta) - 1) * (-teta))) - (-1/teta) * 
        ((1 - p2)^((-teta) - 1) + (1 - p2)^((-teta) - 1) * log((1 - 
            p2)) * (-teta) - (p1^(-teta) * ((1 - p2)^((-teta) - 
            1) + (1 - p2)^((-teta) - 1) * log((1 - p2)) * (-teta)) + 
            p1^(-teta) * log(p1) * ((1 - p2)^((-teta) - 1) * 
                (-teta)))))))*(-exp(teta.st))



bit1.th2 <-  ((p1^((exp(teta.st) + 1 + epsilon)) + (1 - p2)^((exp(teta.st) + 
    1 + epsilon)) - p1^((exp(teta.st) + 1 + epsilon)) * (1 - 
    p2)^((exp(teta.st) + 1 + epsilon)))^(((1/(exp(teta.st) + 
    1 + epsilon)) - 1) - 1) * (((1/(exp(teta.st) + 1 + epsilon)) - 
    1) * (p1^((exp(teta.st) + 1 + epsilon)) * (log(p1) * exp(teta.st)) + 
    (1 - p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - p2)) * 
        exp(teta.st)) - (p1^((exp(teta.st) + 1 + epsilon)) * 
    (log(p1) * exp(teta.st)) * (1 - p2)^((exp(teta.st) + 1 + 
    epsilon)) + p1^((exp(teta.st) + 1 + epsilon)) * ((1 - p2)^((exp(teta.st) + 
    1 + epsilon)) * (log((1 - p2)) * exp(teta.st)))))) - (p1^((exp(teta.st) + 
    1 + epsilon)) + (1 - p2)^((exp(teta.st) + 1 + epsilon)) - 
    p1^((exp(teta.st) + 1 + epsilon)) * (1 - p2)^((exp(teta.st) + 
        1 + epsilon)))^((1/(exp(teta.st) + 1 + epsilon)) - 1) * 
    (log((p1^((exp(teta.st) + 1 + epsilon)) + (1 - p2)^((exp(teta.st) + 
        1 + epsilon)) - p1^((exp(teta.st) + 1 + epsilon)) * (1 - 
        p2)^((exp(teta.st) + 1 + epsilon)))) * (exp(teta.st)/(exp(teta.st) + 
        1 + epsilon)^2))) * ((1/(exp(teta.st) + 1 + epsilon)) * 
    (p1^((exp(teta.st) + 1 + epsilon)) * (log(p1) * exp(teta.st)) + 
        (1 - p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - p2)) * 
            exp(teta.st)) - (p1^((exp(teta.st) + 1 + epsilon)) * 
        (log(p1) * exp(teta.st)) * (1 - p2)^((exp(teta.st) + 
        1 + epsilon)) + p1^((exp(teta.st) + 1 + epsilon)) * ((1 - 
        p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - p2)) * 
        exp(teta.st)))))) + (p1^((exp(teta.st) + 1 + epsilon)) + 
    (1 - p2)^((exp(teta.st) + 1 + epsilon)) - p1^((exp(teta.st) + 
    1 + epsilon)) * (1 - p2)^((exp(teta.st) + 1 + epsilon)))^((1/(exp(teta.st) + 
    1 + epsilon)) - 1) * ((1/(exp(teta.st) + 1 + epsilon)) * 
    (p1^((exp(teta.st) + 1 + epsilon)) * (log(p1) * exp(teta.st)) * 
        (log(p1) * exp(teta.st)) + p1^((exp(teta.st) + 1 + epsilon)) * 
        (log(p1) * exp(teta.st)) + ((1 - p2)^((exp(teta.st) + 
        1 + epsilon)) * (log((1 - p2)) * exp(teta.st)) * (log((1 - 
        p2)) * exp(teta.st)) + (1 - p2)^((exp(teta.st) + 1 + 
        epsilon)) * (log((1 - p2)) * exp(teta.st))) - ((p1^((exp(teta.st) + 
        1 + epsilon)) * (log(p1) * exp(teta.st)) * (log(p1) * 
        exp(teta.st)) + p1^((exp(teta.st) + 1 + epsilon)) * (log(p1) * 
        exp(teta.st))) * (1 - p2)^((exp(teta.st) + 1 + epsilon)) + 
        p1^((exp(teta.st) + 1 + epsilon)) * (log(p1) * exp(teta.st)) * 
            ((1 - p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - 
                p2)) * exp(teta.st))) + (p1^((exp(teta.st) + 
        1 + epsilon)) * (log(p1) * exp(teta.st)) * ((1 - p2)^((exp(teta.st) + 
        1 + epsilon)) * (log((1 - p2)) * exp(teta.st))) + p1^((exp(teta.st) + 
        1 + epsilon)) * ((1 - p2)^((exp(teta.st) + 1 + epsilon)) * 
        (log((1 - p2)) * exp(teta.st)) * (log((1 - p2)) * exp(teta.st)) + 
        (1 - p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - p2)) * 
            exp(teta.st)))))) - exp(teta.st)/(exp(teta.st) + 
    1 + epsilon)^2 * (p1^((exp(teta.st) + 1 + epsilon)) * (log(p1) * 
    exp(teta.st)) + (1 - p2)^((exp(teta.st) + 1 + epsilon)) * 
    (log((1 - p2)) * exp(teta.st)) - (p1^((exp(teta.st) + 1 + 
    epsilon)) * (log(p1) * exp(teta.st)) * (1 - p2)^((exp(teta.st) + 
    1 + epsilon)) + p1^((exp(teta.st) + 1 + epsilon)) * ((1 - 
    p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - p2)) * exp(teta.st)))))) - 
    (((p1^((exp(teta.st) + 1 + epsilon)) + (1 - p2)^((exp(teta.st) + 
        1 + epsilon)) - p1^((exp(teta.st) + 1 + epsilon)) * (1 - 
        p2)^((exp(teta.st) + 1 + epsilon)))^((1/(exp(teta.st) + 
        1 + epsilon)) - 1) * ((1/(exp(teta.st) + 1 + epsilon)) * 
        (p1^((exp(teta.st) + 1 + epsilon)) * (log(p1) * exp(teta.st)) + 
            (1 - p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - 
                p2)) * exp(teta.st)) - (p1^((exp(teta.st) + 1 + 
            epsilon)) * (log(p1) * exp(teta.st)) * (1 - p2)^((exp(teta.st) + 
            1 + epsilon)) + p1^((exp(teta.st) + 1 + epsilon)) * 
            ((1 - p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - 
                p2)) * exp(teta.st)))))) - (p1^((exp(teta.st) + 
        1 + epsilon)) + (1 - p2)^((exp(teta.st) + 1 + epsilon)) - 
        p1^((exp(teta.st) + 1 + epsilon)) * (1 - p2)^((exp(teta.st) + 
            1 + epsilon)))^(1/(exp(teta.st) + 1 + epsilon)) * 
        (log((p1^((exp(teta.st) + 1 + epsilon)) + (1 - p2)^((exp(teta.st) + 
            1 + epsilon)) - p1^((exp(teta.st) + 1 + epsilon)) * 
            (1 - p2)^((exp(teta.st) + 1 + epsilon)))) * (exp(teta.st)/(exp(teta.st) + 
            1 + epsilon)^2))) * (log((p1^((exp(teta.st) + 1 + 
        epsilon)) + (1 - p2)^((exp(teta.st) + 1 + epsilon)) - 
        p1^((exp(teta.st) + 1 + epsilon)) * (1 - p2)^((exp(teta.st) + 
            1 + epsilon)))) * (exp(teta.st)/(exp(teta.st) + 1 + 
        epsilon)^2)) + (p1^((exp(teta.st) + 1 + epsilon)) + (1 - 
        p2)^((exp(teta.st) + 1 + epsilon)) - p1^((exp(teta.st) + 
        1 + epsilon)) * (1 - p2)^((exp(teta.st) + 1 + epsilon)))^(1/(exp(teta.st) + 
        1 + epsilon)) * ((p1^((exp(teta.st) + 1 + epsilon)) * 
        (log(p1) * exp(teta.st)) + (1 - p2)^((exp(teta.st) + 
        1 + epsilon)) * (log((1 - p2)) * exp(teta.st)) - (p1^((exp(teta.st) + 
        1 + epsilon)) * (log(p1) * exp(teta.st)) * (1 - p2)^((exp(teta.st) + 
        1 + epsilon)) + p1^((exp(teta.st) + 1 + epsilon)) * ((1 - 
        p2)^((exp(teta.st) + 1 + epsilon)) * (log((1 - p2)) * 
        exp(teta.st)))))/(p1^((exp(teta.st) + 1 + epsilon)) + 
        (1 - p2)^((exp(teta.st) + 1 + epsilon)) - p1^((exp(teta.st) + 
        1 + epsilon)) * (1 - p2)^((exp(teta.st) + 1 + epsilon))) * 
        (exp(teta.st)/(exp(teta.st) + 1 + epsilon)^2) + log((p1^((exp(teta.st) + 
        1 + epsilon)) + (1 - p2)^((exp(teta.st) + 1 + epsilon)) - 
        p1^((exp(teta.st) + 1 + epsilon)) * (1 - p2)^((exp(teta.st) + 
            1 + epsilon)))) * (exp(teta.st)/(exp(teta.st) + 1 + 
        epsilon)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
        1 + epsilon)))/((exp(teta.st) + 1 + epsilon)^2)^2)))





}


if(BivD=="J180"){


  
 c.copula.be1 <- 1 - (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 1) * 
     ((1/teta) * (p1^(teta - 1) * teta - p1^(teta - 1) * teta * 
         p2^teta))
 
 
 
  
   c.copula.be2 <- 1 - (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 1) * 
     ((1/teta) * (p2^(teta - 1) * teta - p1^teta * (p2^(teta - 
         1) * teta)))
 
 
 
 
   c.copula.theta <-   (-((p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 1) * ((1/teta) * 
     (p1^teta * log(p1) + p2^teta * log(p2) - (p1^teta * log(p1) * 
         p2^teta + p1^teta * (p2^teta * log(p2))))) - (p1^teta + 
     p2^teta - p1^teta * p2^teta)^(1/teta) * (log((p1^teta + p2^teta - 
    p1^teta * p2^teta)) * (1/teta^2))))*exp(teta.st)
    
c.copula2.be1 <- -((p1^teta + p2^teta - p1^teta * p2^teta)^(((1/teta) - 1) - 1) * 
    (((1/teta) - 1) * (p1^(teta - 1) * teta - p1^(teta - 1) * 
        teta * p2^teta)) * ((1/teta) * (p1^(teta - 1) * teta - 
    p1^(teta - 1) * teta * p2^teta)) + (p1^teta + p2^teta - p1^teta * 
    p2^teta)^((1/teta) - 1) * ((1/teta) * (p1^((teta - 1) - 1) * 
    (teta - 1) * teta - p1^((teta - 1) - 1) * (teta - 1) * teta * 
    p2^teta)))




                  
c.copula2.be2 <- -((p1^teta + p2^teta - p1^teta * p2^teta)^(((1/teta) - 1) - 1) * 
    (((1/teta) - 1) * (p2^(teta - 1) * teta - p1^teta * (p2^(teta - 
        1) * teta))) * ((1/teta) * (p2^(teta - 1) * teta - p1^teta * 
    (p2^(teta - 1) * teta))) + (p1^teta + p2^teta - p1^teta * 
    p2^teta)^((1/teta) - 1) * ((1/teta) * (p2^((teta - 1) - 1) * 
    (teta - 1) * teta - p1^teta * (p2^((teta - 1) - 1) * (teta - 
    1) * teta))))



c.copula2.be1be2 <- -((p1^teta + p2^teta - p1^teta * p2^teta)^(((1/teta) - 1) - 1) * 
    (((1/teta) - 1) * (p2^(teta - 1) * teta - p1^teta * (p2^(teta - 
        1) * teta))) * ((1/teta) * (p1^(teta - 1) * teta - p1^(teta - 
    1) * teta * p2^teta)) - (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 
    1) * ((1/teta) * (p1^(teta - 1) * teta * (p2^(teta - 1) * 
    teta))))

c.copula2.be1th <-(-(((p1^teta + p2^teta - p1^teta * p2^teta)^(((1/teta) - 1) - 
    1) * (((1/teta) - 1) * (p1^teta * log(p1) + p2^teta * log(p2) - 
    (p1^teta * log(p1) * p2^teta + p1^teta * (p2^teta * log(p2))))) - 
    (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 1) * 
        (log((p1^teta + p2^teta - p1^teta * p2^teta)) * (1/teta^2))) * 
    ((1/teta) * (p1^(teta - 1) * teta - p1^(teta - 1) * teta * 
        p2^teta)) + (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 
    1) * ((1/teta) * (p1^(teta - 1) * log(p1) * teta + p1^(teta - 
    1) - ((p1^(teta - 1) * log(p1) * teta + p1^(teta - 1)) * 
    p2^teta + p1^(teta - 1) * teta * (p2^teta * log(p2)))) - 
    1/teta^2 * (p1^(teta - 1) * teta - p1^(teta - 1) * teta * 
        p2^teta))))*exp(teta.st)

c.copula2.be2th <- (-(((p1^teta + p2^teta - p1^teta * p2^teta)^(((1/teta) - 1) - 
    1) * (((1/teta) - 1) * (p1^teta * log(p1) + p2^teta * log(p2) - 
    (p1^teta * log(p1) * p2^teta + p1^teta * (p2^teta * log(p2))))) - 
    (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 1) * 
        (log((p1^teta + p2^teta - p1^teta * p2^teta)) * (1/teta^2))) * 
    ((1/teta) * (p2^(teta - 1) * teta - p1^teta * (p2^(teta - 
        1) * teta))) + (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 
    1) * ((1/teta) * (p2^(teta - 1) * log(p2) * teta + p2^(teta - 
    1) - (p1^teta * log(p1) * (p2^(teta - 1) * teta) + p1^teta * 
    (p2^(teta - 1) * log(p2) * teta + p2^(teta - 1)))) - 1/teta^2 * 
    (p2^(teta - 1) * teta - p1^teta * (p2^(teta - 1) * teta)))))*exp(teta.st)



bit1.th2 <- -(((p1^(exp(teta.st) + 1 + epsilon) + p2^(exp(teta.st) + 1 + 
    epsilon) - p1^(exp(teta.st) + 1 + epsilon) * p2^(exp(teta.st) + 
    1 + epsilon))^(((1/(exp(teta.st) + 1 + epsilon)) - 1) - 1) * 
    (((1/(exp(teta.st) + 1 + epsilon)) - 1) * (p1^(exp(teta.st) + 
        1 + epsilon) * (log(p1) * exp(teta.st)) + p2^(exp(teta.st) + 
        1 + epsilon) * (log(p2) * exp(teta.st)) - (p1^(exp(teta.st) + 
        1 + epsilon) * (log(p1) * exp(teta.st)) * p2^(exp(teta.st) + 
        1 + epsilon) + p1^(exp(teta.st) + 1 + epsilon) * (p2^(exp(teta.st) + 
        1 + epsilon) * (log(p2) * exp(teta.st)))))) - (p1^(exp(teta.st) + 
    1 + epsilon) + p2^(exp(teta.st) + 1 + epsilon) - p1^(exp(teta.st) + 
    1 + epsilon) * p2^(exp(teta.st) + 1 + epsilon))^((1/(exp(teta.st) + 
    1 + epsilon)) - 1) * (log((p1^(exp(teta.st) + 1 + epsilon) + 
    p2^(exp(teta.st) + 1 + epsilon) - p1^(exp(teta.st) + 1 + 
    epsilon) * p2^(exp(teta.st) + 1 + epsilon))) * (exp(teta.st)/(exp(teta.st) + 
    1 + epsilon)^2))) * ((1/(exp(teta.st) + 1 + epsilon)) * (p1^(exp(teta.st) + 
    1 + epsilon) * (log(p1) * exp(teta.st)) + p2^(exp(teta.st) + 
    1 + epsilon) * (log(p2) * exp(teta.st)) - (p1^(exp(teta.st) + 
    1 + epsilon) * (log(p1) * exp(teta.st)) * p2^(exp(teta.st) + 
    1 + epsilon) + p1^(exp(teta.st) + 1 + epsilon) * (p2^(exp(teta.st) + 
    1 + epsilon) * (log(p2) * exp(teta.st)))))) + (p1^(exp(teta.st) + 
    1 + epsilon) + p2^(exp(teta.st) + 1 + epsilon) - p1^(exp(teta.st) + 
    1 + epsilon) * p2^(exp(teta.st) + 1 + epsilon))^((1/(exp(teta.st) + 
    1 + epsilon)) - 1) * ((1/(exp(teta.st) + 1 + epsilon)) * 
    (p1^(exp(teta.st) + 1 + epsilon) * (log(p1) * exp(teta.st)) * 
        (log(p1) * exp(teta.st)) + p1^(exp(teta.st) + 1 + epsilon) * 
        (log(p1) * exp(teta.st)) + (p2^(exp(teta.st) + 1 + epsilon) * 
        (log(p2) * exp(teta.st)) * (log(p2) * exp(teta.st)) + 
        p2^(exp(teta.st) + 1 + epsilon) * (log(p2) * exp(teta.st))) - 
        ((p1^(exp(teta.st) + 1 + epsilon) * (log(p1) * exp(teta.st)) * 
            (log(p1) * exp(teta.st)) + p1^(exp(teta.st) + 1 + 
            epsilon) * (log(p1) * exp(teta.st))) * p2^(exp(teta.st) + 
            1 + epsilon) + p1^(exp(teta.st) + 1 + epsilon) * 
            (log(p1) * exp(teta.st)) * (p2^(exp(teta.st) + 1 + 
            epsilon) * (log(p2) * exp(teta.st))) + (p1^(exp(teta.st) + 
            1 + epsilon) * (log(p1) * exp(teta.st)) * (p2^(exp(teta.st) + 
            1 + epsilon) * (log(p2) * exp(teta.st))) + p1^(exp(teta.st) + 
            1 + epsilon) * (p2^(exp(teta.st) + 1 + epsilon) * 
            (log(p2) * exp(teta.st)) * (log(p2) * exp(teta.st)) + 
            p2^(exp(teta.st) + 1 + epsilon) * (log(p2) * exp(teta.st)))))) - 
    exp(teta.st)/(exp(teta.st) + 1 + epsilon)^2 * (p1^(exp(teta.st) + 
        1 + epsilon) * (log(p1) * exp(teta.st)) + p2^(exp(teta.st) + 
        1 + epsilon) * (log(p2) * exp(teta.st)) - (p1^(exp(teta.st) + 
        1 + epsilon) * (log(p1) * exp(teta.st)) * p2^(exp(teta.st) + 
        1 + epsilon) + p1^(exp(teta.st) + 1 + epsilon) * (p2^(exp(teta.st) + 
        1 + epsilon) * (log(p2) * exp(teta.st)))))) - (((p1^(exp(teta.st) + 
    1 + epsilon) + p2^(exp(teta.st) + 1 + epsilon) - p1^(exp(teta.st) + 
    1 + epsilon) * p2^(exp(teta.st) + 1 + epsilon))^((1/(exp(teta.st) + 
    1 + epsilon)) - 1) * ((1/(exp(teta.st) + 1 + epsilon)) * 
    (p1^(exp(teta.st) + 1 + epsilon) * (log(p1) * exp(teta.st)) + 
        p2^(exp(teta.st) + 1 + epsilon) * (log(p2) * exp(teta.st)) - 
        (p1^(exp(teta.st) + 1 + epsilon) * (log(p1) * exp(teta.st)) * 
            p2^(exp(teta.st) + 1 + epsilon) + p1^(exp(teta.st) + 
            1 + epsilon) * (p2^(exp(teta.st) + 1 + epsilon) * 
            (log(p2) * exp(teta.st)))))) - (p1^(exp(teta.st) + 
    1 + epsilon) + p2^(exp(teta.st) + 1 + epsilon) - p1^(exp(teta.st) + 
    1 + epsilon) * p2^(exp(teta.st) + 1 + epsilon))^(1/(exp(teta.st) + 
    1 + epsilon)) * (log((p1^(exp(teta.st) + 1 + epsilon) + p2^(exp(teta.st) + 
    1 + epsilon) - p1^(exp(teta.st) + 1 + epsilon) * p2^(exp(teta.st) + 
    1 + epsilon))) * (exp(teta.st)/(exp(teta.st) + 1 + epsilon)^2))) * 
    (log((p1^(exp(teta.st) + 1 + epsilon) + p2^(exp(teta.st) + 
        1 + epsilon) - p1^(exp(teta.st) + 1 + epsilon) * p2^(exp(teta.st) + 
        1 + epsilon))) * (exp(teta.st)/(exp(teta.st) + 1 + epsilon)^2)) + 
    (p1^(exp(teta.st) + 1 + epsilon) + p2^(exp(teta.st) + 1 + 
        epsilon) - p1^(exp(teta.st) + 1 + epsilon) * p2^(exp(teta.st) + 
        1 + epsilon))^(1/(exp(teta.st) + 1 + epsilon)) * ((p1^(exp(teta.st) + 
        1 + epsilon) * (log(p1) * exp(teta.st)) + p2^(exp(teta.st) + 
        1 + epsilon) * (log(p2) * exp(teta.st)) - (p1^(exp(teta.st) + 
        1 + epsilon) * (log(p1) * exp(teta.st)) * p2^(exp(teta.st) + 
        1 + epsilon) + p1^(exp(teta.st) + 1 + epsilon) * (p2^(exp(teta.st) + 
        1 + epsilon) * (log(p2) * exp(teta.st)))))/(p1^(exp(teta.st) + 
        1 + epsilon) + p2^(exp(teta.st) + 1 + epsilon) - p1^(exp(teta.st) + 
        1 + epsilon) * p2^(exp(teta.st) + 1 + epsilon)) * (exp(teta.st)/(exp(teta.st) + 
        1 + epsilon)^2) + log((p1^(exp(teta.st) + 1 + epsilon) + 
        p2^(exp(teta.st) + 1 + epsilon) - p1^(exp(teta.st) + 
        1 + epsilon) * p2^(exp(teta.st) + 1 + epsilon))) * (exp(teta.st)/(exp(teta.st) + 
        1 + epsilon)^2 - exp(teta.st) * (2 * (exp(teta.st) * 
        (exp(teta.st) + 1 + epsilon)))/((exp(teta.st) + 1 + epsilon)^2)^2))))




}



if(BivD=="J270"){



  c.copula.be1 <- 1 - ((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - p1)^((-teta) - 1) * (-teta) - (1 - 
    p1)^((-teta) - 1) * (-teta) * p2^(-teta)))


 
  c.copula.be2 <- ((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta) - (1 - p1)^(-teta) * 
    (p2^((-teta) - 1) * (-teta))))


  c.copula.theta <-  (((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^(-1/teta) * 
    (log(((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * 
        p2^(-teta))) * (1/teta^2)) - ((1 - p1)^(-teta) + p2^(-teta) - 
    (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 1) * ((-1/teta) * 
    (p2^(-teta) * log(p2) + (1 - p1)^(-teta) * log((1 - p1)) - 
        ((1 - p1)^(-teta) * (p2^(-teta) * log(p2)) + (1 - p1)^(-teta) * 
            log((1 - p1)) * p2^(-teta)))))*(-exp(teta.st))
            
c.copula2.be1 <- ((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - p1)^(((-teta) - 1) - 1) * ((-teta) - 
    1) * (-teta) - (1 - p1)^(((-teta) - 1) - 1) * ((-teta) - 
    1) * (-teta) * p2^(-teta))) + ((1 - p1)^(-teta) + p2^(-teta) - 
    (1 - p1)^(-teta) * p2^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * ((1 - p1)^((-teta) - 1) * (-teta) - (1 - p1)^((-teta) - 
    1) * (-teta) * p2^(-teta))) * ((-1/teta) * ((1 - p1)^((-teta) - 
    1) * (-teta) - (1 - p1)^((-teta) - 1) * (-teta) * p2^(-teta)))



                  
c.copula2.be2 <- ((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (p2^((-teta) - 1) * (-teta) - 
    (1 - p1)^(-teta) * (p2^((-teta) - 1) * (-teta)))) * ((-1/teta) * 
    (p2^((-teta) - 1) * (-teta) - (1 - p1)^(-teta) * (p2^((-teta) - 
        1) * (-teta)))) + ((1 - p1)^(-teta) + p2^(-teta) - (1 - 
    p1)^(-teta) * p2^(-teta))^((-1/teta) - 1) * ((-1/teta) * 
    (p2^(((-teta) - 1) - 1) * ((-teta) - 1) * (-teta) - (1 - 
        p1)^(-teta) * (p2^(((-teta) - 1) - 1) * ((-teta) - 1) * 
        (-teta))))


c.copula2.be1be2 <- -(((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (p2^((-teta) - 1) * (-teta) - 
    (1 - p1)^(-teta) * (p2^((-teta) - 1) * (-teta)))) * ((-1/teta) * 
    ((1 - p1)^((-teta) - 1) * (-teta) - (1 - p1)^((-teta) - 1) * 
        (-teta) * p2^(-teta))) - ((1 - p1)^(-teta) + p2^(-teta) - 
    (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 1) * ((-1/teta) * 
    ((1 - p1)^((-teta) - 1) * (-teta) * (p2^((-teta) - 1) * (-teta)))))

c.copula2.be1th <-(-((((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 
    1) * (log(((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * 
    p2^(-teta))) * (1/teta^2)) - ((1 - p1)^(-teta) + p2^(-teta) - 
    (1 - p1)^(-teta) * p2^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (p2^(-teta) * log(p2) + (1 - p1)^(-teta) * log((1 - 
    p1)) - ((1 - p1)^(-teta) * (p2^(-teta) * log(p2)) + (1 - 
    p1)^(-teta) * log((1 - p1)) * p2^(-teta))))) * ((-1/teta) * 
    ((1 - p1)^((-teta) - 1) * (-teta) - (1 - p1)^((-teta) - 1) * 
        (-teta) * p2^(-teta))) + ((1 - p1)^(-teta) + p2^(-teta) - 
    (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 1) * (1/teta^2 * 
    ((1 - p1)^((-teta) - 1) * (-teta) - (1 - p1)^((-teta) - 1) * 
        (-teta) * p2^(-teta)) - (-1/teta) * ((1 - p1)^((-teta) - 
    1) + (1 - p1)^((-teta) - 1) * log((1 - p1)) * (-teta) - ((1 - 
    p1)^((-teta) - 1) * (-teta) * (p2^(-teta) * log(p2)) + ((1 - 
    p1)^((-teta) - 1) + (1 - p1)^((-teta) - 1) * log((1 - p1)) * 
    (-teta)) * p2^(-teta))))))*(-exp(teta.st))

c.copula2.be2th <- ((((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 
    1) * (log(((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * 
    p2^(-teta))) * (1/teta^2)) - ((1 - p1)^(-teta) + p2^(-teta) - 
    (1 - p1)^(-teta) * p2^(-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (p2^(-teta) * log(p2) + (1 - p1)^(-teta) * log((1 - 
    p1)) - ((1 - p1)^(-teta) * (p2^(-teta) * log(p2)) + (1 - 
    p1)^(-teta) * log((1 - p1)) * p2^(-teta))))) * ((-1/teta) * 
    (p2^((-teta) - 1) * (-teta) - (1 - p1)^(-teta) * (p2^((-teta) - 
        1) * (-teta)))) + ((1 - p1)^(-teta) + p2^(-teta) - (1 - 
    p1)^(-teta) * p2^(-teta))^((-1/teta) - 1) * (1/teta^2 * (p2^((-teta) - 
    1) * (-teta) - (1 - p1)^(-teta) * (p2^((-teta) - 1) * (-teta))) - 
    (-1/teta) * (p2^((-teta) - 1) + p2^((-teta) - 1) * log(p2) * 
        (-teta) - ((1 - p1)^(-teta) * (p2^((-teta) - 1) + p2^((-teta) - 
        1) * log(p2) * (-teta)) + (1 - p1)^(-teta) * log((1 - 
        p1)) * (p2^((-teta) - 1) * (-teta))))))*(-exp(teta.st))




bit1.th2 <- (((1 - p1)^((exp(teta.st) + 1 + epsilon)) + p2^((exp(teta.st) + 
    1 + epsilon)) - (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
    p2^((exp(teta.st) + 1 + epsilon)))^(((1/((exp(teta.st) + 
    1 + epsilon))) - 1) - 1) * (((1/((exp(teta.st) + 1 + epsilon))) - 
    1) * ((1 - p1)^((exp(teta.st) + 1 + epsilon)) * (log((1 - 
    p1)) * exp(teta.st)) + p2^((exp(teta.st) + 1 + epsilon)) * 
    (log(p2) * exp(teta.st)) - ((1 - p1)^((exp(teta.st) + 1 + 
    epsilon)) * (log((1 - p1)) * exp(teta.st)) * p2^((exp(teta.st) + 
    1 + epsilon)) + (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
    (p2^((exp(teta.st) + 1 + epsilon)) * (log(p2) * exp(teta.st)))))) - 
    ((1 - p1)^((exp(teta.st) + 1 + epsilon)) + p2^((exp(teta.st) + 
        1 + epsilon)) - (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        p2^((exp(teta.st) + 1 + epsilon)))^((1/((exp(teta.st) + 
        1 + epsilon))) - 1) * (log(((1 - p1)^((exp(teta.st) + 
        1 + epsilon)) + p2^((exp(teta.st) + 1 + epsilon)) - (1 - 
        p1)^((exp(teta.st) + 1 + epsilon)) * p2^((exp(teta.st) + 
        1 + epsilon)))) * (exp(teta.st)/((exp(teta.st) + 1 + 
        epsilon))^2))) * ((1/((exp(teta.st) + 1 + epsilon))) * 
    ((1 - p1)^((exp(teta.st) + 1 + epsilon)) * (log((1 - p1)) * 
        exp(teta.st)) + p2^((exp(teta.st) + 1 + epsilon)) * (log(p2) * 
        exp(teta.st)) - ((1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        (log((1 - p1)) * exp(teta.st)) * p2^((exp(teta.st) + 
        1 + epsilon)) + (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        (p2^((exp(teta.st) + 1 + epsilon)) * (log(p2) * exp(teta.st)))))) + 
    ((1 - p1)^((exp(teta.st) + 1 + epsilon)) + p2^((exp(teta.st) + 
        1 + epsilon)) - (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        p2^((exp(teta.st) + 1 + epsilon)))^((1/((exp(teta.st) + 
        1 + epsilon))) - 1) * ((1/((exp(teta.st) + 1 + epsilon))) * 
        ((1 - p1)^((exp(teta.st) + 1 + epsilon)) * (log((1 - 
            p1)) * exp(teta.st)) * (log((1 - p1)) * exp(teta.st)) + 
            (1 - p1)^((exp(teta.st) + 1 + epsilon)) * (log((1 - 
                p1)) * exp(teta.st)) + (p2^((exp(teta.st) + 1 + 
            epsilon)) * (log(p2) * exp(teta.st)) * (log(p2) * 
            exp(teta.st)) + p2^((exp(teta.st) + 1 + epsilon)) * 
            (log(p2) * exp(teta.st))) - (((1 - p1)^((exp(teta.st) + 
            1 + epsilon)) * (log((1 - p1)) * exp(teta.st)) * 
            (log((1 - p1)) * exp(teta.st)) + (1 - p1)^((exp(teta.st) + 
            1 + epsilon)) * (log((1 - p1)) * exp(teta.st))) * 
            p2^((exp(teta.st) + 1 + epsilon)) + (1 - p1)^((exp(teta.st) + 
            1 + epsilon)) * (log((1 - p1)) * exp(teta.st)) * 
            (p2^((exp(teta.st) + 1 + epsilon)) * (log(p2) * exp(teta.st))) + 
            ((1 - p1)^((exp(teta.st) + 1 + epsilon)) * (log((1 - 
                p1)) * exp(teta.st)) * (p2^((exp(teta.st) + 1 + 
                epsilon)) * (log(p2) * exp(teta.st))) + (1 - 
                p1)^((exp(teta.st) + 1 + epsilon)) * (p2^((exp(teta.st) + 
                1 + epsilon)) * (log(p2) * exp(teta.st)) * (log(p2) * 
                exp(teta.st)) + p2^((exp(teta.st) + 1 + epsilon)) * 
                (log(p2) * exp(teta.st)))))) - exp(teta.st)/((exp(teta.st) + 
        1 + epsilon))^2 * ((1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        (log((1 - p1)) * exp(teta.st)) + p2^((exp(teta.st) + 
        1 + epsilon)) * (log(p2) * exp(teta.st)) - ((1 - p1)^((exp(teta.st) + 
        1 + epsilon)) * (log((1 - p1)) * exp(teta.st)) * p2^((exp(teta.st) + 
        1 + epsilon)) + (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        (p2^((exp(teta.st) + 1 + epsilon)) * (log(p2) * exp(teta.st)))))) - 
    ((((1 - p1)^((exp(teta.st) + 1 + epsilon)) + p2^((exp(teta.st) + 
        1 + epsilon)) - (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        p2^((exp(teta.st) + 1 + epsilon)))^((1/((exp(teta.st) + 
        1 + epsilon))) - 1) * ((1/((exp(teta.st) + 1 + epsilon))) * 
        ((1 - p1)^((exp(teta.st) + 1 + epsilon)) * (log((1 - 
            p1)) * exp(teta.st)) + p2^((exp(teta.st) + 1 + epsilon)) * 
            (log(p2) * exp(teta.st)) - ((1 - p1)^((exp(teta.st) + 
            1 + epsilon)) * (log((1 - p1)) * exp(teta.st)) * 
            p2^((exp(teta.st) + 1 + epsilon)) + (1 - p1)^((exp(teta.st) + 
            1 + epsilon)) * (p2^((exp(teta.st) + 1 + epsilon)) * 
            (log(p2) * exp(teta.st)))))) - ((1 - p1)^((exp(teta.st) + 
        1 + epsilon)) + p2^((exp(teta.st) + 1 + epsilon)) - (1 - 
        p1)^((exp(teta.st) + 1 + epsilon)) * p2^((exp(teta.st) + 
        1 + epsilon)))^(1/((exp(teta.st) + 1 + epsilon))) * (log(((1 - 
        p1)^((exp(teta.st) + 1 + epsilon)) + p2^((exp(teta.st) + 
        1 + epsilon)) - (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        p2^((exp(teta.st) + 1 + epsilon)))) * (exp(teta.st)/((exp(teta.st) + 
        1 + epsilon))^2))) * (log(((1 - p1)^((exp(teta.st) + 
        1 + epsilon)) + p2^((exp(teta.st) + 1 + epsilon)) - (1 - 
        p1)^((exp(teta.st) + 1 + epsilon)) * p2^((exp(teta.st) + 
        1 + epsilon)))) * (exp(teta.st)/((exp(teta.st) + 1 + 
        epsilon))^2)) + ((1 - p1)^((exp(teta.st) + 1 + epsilon)) + 
        p2^((exp(teta.st) + 1 + epsilon)) - (1 - p1)^((exp(teta.st) + 
        1 + epsilon)) * p2^((exp(teta.st) + 1 + epsilon)))^(1/((exp(teta.st) + 
        1 + epsilon))) * (((1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        (log((1 - p1)) * exp(teta.st)) + p2^((exp(teta.st) + 
        1 + epsilon)) * (log(p2) * exp(teta.st)) - ((1 - p1)^((exp(teta.st) + 
        1 + epsilon)) * (log((1 - p1)) * exp(teta.st)) * p2^((exp(teta.st) + 
        1 + epsilon)) + (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        (p2^((exp(teta.st) + 1 + epsilon)) * (log(p2) * exp(teta.st)))))/((1 - 
        p1)^((exp(teta.st) + 1 + epsilon)) + p2^((exp(teta.st) + 
        1 + epsilon)) - (1 - p1)^((exp(teta.st) + 1 + epsilon)) * 
        p2^((exp(teta.st) + 1 + epsilon))) * (exp(teta.st)/((exp(teta.st) + 
        1 + epsilon))^2) + log(((1 - p1)^((exp(teta.st) + 1 + 
        epsilon)) + p2^((exp(teta.st) + 1 + epsilon)) - (1 - 
        p1)^((exp(teta.st) + 1 + epsilon)) * p2^((exp(teta.st) + 
        1 + epsilon)))) * (exp(teta.st)/((exp(teta.st) + 1 + 
        epsilon))^2 - exp(teta.st) * (2 * (exp(teta.st) * ((exp(teta.st) + 
        1 + epsilon))))/(((exp(teta.st) + 1 + epsilon))^2)^2)))




}







         list(c.copula.be1=c.copula.be1, c.copula.be2=c.copula.be2, c.copula.theta=c.copula.theta, 
              c.copula2.be1=c.copula2.be1, c.copula2.be2=c.copula2.be2, c.copula2.be1be2=c.copula2.be1be2,
              c.copula2.be1th=c.copula2.be1th,c.copula2.be2th=c.copula2.be2th,bit1.th2=bit1.th2)     


 

}




     
























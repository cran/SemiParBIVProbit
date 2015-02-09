copgHs <- function(p1,p2,eta1=NULL,eta2=NULL,teta,teta.st,xi1,xi1.st=NULL,xi2,xi2.st=NULL,BivD,nC,nu=NULL,PL=NULL,eqPL=NULL){

epsilon <- 0 # .Machine$double.eps*10^6
c.copula.lambda1 <- c.copula.lambda2 <- bit1.lambda1.2 <- bit1.lambda2.2 <- c.copula2.be1lambda1 <- c.copula2.be2lambda2 <- c.copula2.be1lambda2 <- c.copula2.be2lambda1 <- bit1.thlambda1 <- bit1.thlambda2 <- bit1.lambda1lambda2 <- NULL
der.p1.lambda1 <- der.p2.lambda2 <- der.d.n1.be1 <- der.d.n2.be2 <- der.der.p1.lam1.der.p1 <- der.der.p2.lam2.der.p2 <- der2.p1.lambda1 <- der2.p2.lambda2 <- der.d.n1.lambda1 <- der.d.n2.lambda2 <- der.p1.lambda1 <- der.p2.lambda2 <- NULL
eps <- 1e-4

if(BivD=="N"){

tt.st <- tanh(teta.st)
qp1 <- qnorm(p1)
qp2 <- qnorm(p2)
ct.st <- cosh(teta.st)

c.copula.be1 <- pnorm( (qp2 - teta*qp1)/sqrt(1 - teta^2)   )  # BiCopHfunc(p1, p2, family=1, par=teta)$hfunc1                          
c.copula.be2 <- pnorm( (qp1 - teta*qp2)/sqrt(1 - teta^2)   )  # BiCopHfunc(p1, p2, family=1, par=teta)$hfunc2

c.copula.theta <- dbinorm(qp1,qp2, cov12=teta)*(1/ct.st^2) 

c.copula2.be1 <- dnorm((qp2-teta*qp1)/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qp1^2/2) # BiCopHfuncDeriv(p2, p1, 1, par=teta, deriv="u2")    
c.copula2.be2 <- dnorm((qp1-teta*qp2)/sqrt(1 - teta^2))  * sqrt(2*pi) *(-teta)/sqrt(1 - teta^2)/exp(-qp2^2/2) # BiCopHfuncDeriv(p1, p2, 1, par=teta, deriv="u2")

c.copula2.be1be2 <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) # BiCopPDF(p1, p2, 1, par=teta)

c.copula2.be1th <- -(dnorm((qp2 - tt.st * qp1)/sqrt(1 - tt.st^2)) * 
     (1/ct.st^2 * qp1/sqrt(1 - tt.st^2) - (qp2 - 
         tt.st * qp1) * (0.5 * (2 * (1/ct.st^2 * 
         tt.st) * (1 - tt.st^2)^-0.5))/sqrt(1 - 
         tt.st^2)^2)) 

c.copula2.be2th <- -(dnorm((qp1 - tt.st * qp2)/sqrt(1 - tt.st^2)) * 
    (1/ct.st^2 * qp2/sqrt(1 - tt.st^2) - (qp1 - 
        tt.st * qp2) * (0.5 * (2 * (1/ct.st^2 * 
        tt.st) * (1 - tt.st^2)^-0.5))/sqrt(1 - 
        tt.st^2)^2))
 
bit1.th2 <- (2 * pi * (0.5 * (2 * (1/ct.st^2 * tt.st) * (1 - 
    tt.st^2)^-0.5))/(2 * pi * sqrt(1 - tt.st^2))^2 * 
    exp(-1/(2 * (1 - tt.st^2)) * (qp1^2 + qp2^2 - 2 * 
        tt.st * qp1 * qp2)) - 1/(2 * pi * sqrt(1 - 
    tt.st^2)) * (exp(-1/(2 * (1 - tt.st^2)) * 
    (qp1^2 + qp2^2 - 2 * tt.st * qp1 * qp2)) * (-1/(2 * 
    (1 - tt.st^2)) * (2 * (1/ct.st^2) * qp1 * 
    qp2) + 2 * (2 * (1/ct.st^2 * tt.st))/(2 * 
    (1 - tt.st^2))^2 * (qp1^2 + qp2^2 - 2 * tt.st * 
    qp1 * qp2))))/ct.st^2 - 1/(2 * pi * sqrt(1 - tt.st^2)) * 
    exp(-1/(2 * (1 - tt.st^2)) * (qp1^2 + qp2^2 - 2 * 
        tt.st * qp1 * qp2)) * 1 * (2 * (sinh(teta.st) * 
    ct.st))/(ct.st^2)^2 


bit1.th2ATE <- 0.5 * (pi * (0.5 * (2 * teta * (1 - teta^2)^-0.5)))/(pi * sqrt(1 - 
    teta^2))^2 * (exp(-0.5/(1 - teta^2) * (qp1^2 + qp2^2 - 2 * 
    teta * qp1 * qp2))) - 0.5/(pi * sqrt(1 - teta^2)) * (exp(-0.5/(1 - 
    teta^2) * (qp1^2 + qp2^2 - 2 * teta * qp1 * qp2)) * (-0.5/(1 - 
    teta^2) * (2 * qp1 * qp2) + 0.5 * (2 * teta)/(1 - teta^2)^2 * 
    (qp1^2 + qp2^2 - 2 * teta * qp1 * qp2)))




}


if(BivD=="C0"){

et.st <- exp(teta.st)
lo.p1 <- log(p1)
lo.p2 <- log(p2)


c.copula.be1 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p1^((-teta) - 1) * (-teta)))
  c.copula.be2 <- (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * (p2^((-teta) - 1) * (-teta)))
  c.copula.theta <- ((p1^(-teta) + p2^(-teta) - 1)^(-1/teta) * (log((p1^(-teta) + 
    p2^(-teta) - 1)) * (1/teta^2)) - (p1^(-teta) + p2^(-teta) - 
    1)^((-1/teta) - 1) * ((-1/teta) * (p2^(-teta) * lo.p2 + 
    p1^(-teta) * lo.p1)))*et.st

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

c.copula2.be1th <- ((p1^(-(et.st + epsilon)) + p2^(-(et.st + epsilon)) - 
    1)^((-1/(et.st + epsilon)) - 1) * (log((p1^(-(et.st + 
    epsilon)) + p2^(-(et.st + epsilon)) - 1)) * (et.st/(et.st + 
    epsilon)^2)) - (p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^(((-1/(et.st + epsilon)) - 1) - 1) * 
    (((-1/(et.st + epsilon)) - 1) * (p2^(-(et.st + 
        epsilon)) * (lo.p2 * et.st) + p1^(-(et.st + 
        epsilon)) * (lo.p1 * et.st)))) * ((-1/(et.st + 
    epsilon)) * (p1^((-(et.st + epsilon)) - 1) * (-(et.st + 
    epsilon)))) + (p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^((-1/(et.st + epsilon)) - 1) * (et.st/(et.st + 
    epsilon)^2 * (p1^((-(et.st + epsilon)) - 1) * (-(et.st + 
    epsilon))) - (-1/(et.st + epsilon)) * (p1^((-(et.st + 
    epsilon)) - 1) * et.st + p1^((-(et.st + epsilon)) - 
    1) * (lo.p1 * et.st) * (-(et.st + epsilon))))


c.copula2.be2th <- ((p1^(-(et.st + epsilon)) + p2^(-(et.st + epsilon)) - 
    1)^((-1/(et.st + epsilon)) - 1) * (log((p1^(-(et.st + 
    epsilon)) + p2^(-(et.st + epsilon)) - 1)) * (et.st/(et.st + 
    epsilon)^2)) - (p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^(((-1/(et.st + epsilon)) - 1) - 1) * 
    (((-1/(et.st + epsilon)) - 1) * (p2^(-(et.st + 
        epsilon)) * (lo.p2 * et.st) + p1^(-(et.st + 
        epsilon)) * (lo.p1 * et.st)))) * ((-1/(et.st + 
    epsilon)) * (p2^((-(et.st + epsilon)) - 1) * (-(et.st + 
    epsilon)))) + (p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^((-1/(et.st + epsilon)) - 1) * (et.st/(et.st + 
    epsilon)^2 * (p2^((-(et.st + epsilon)) - 1) * (-(et.st + 
    epsilon))) - (-1/(et.st + epsilon)) * (p2^((-(et.st + 
    epsilon)) - 1) * et.st + p2^((-(et.st + epsilon)) - 
    1) * (lo.p2 * et.st) * (-(et.st + epsilon))))

bit1.th2 <-((p1^(-(et.st + epsilon)) + p2^(-(et.st + epsilon)) - 
    1)^(-1/(et.st + epsilon)) * (log((p1^(-(et.st + 
    epsilon)) + p2^(-(et.st + epsilon)) - 1)) * (et.st/(et.st + 
    epsilon)^2)) - (p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^((-1/(et.st + epsilon)) - 1) * ((-1/(et.st + 
    epsilon)) * (p2^(-(et.st + epsilon)) * (lo.p2 * 
    et.st) + p1^(-(et.st + epsilon)) * (lo.p1 * 
    et.st)))) * (log((p1^(-(et.st + epsilon)) + 
    p2^(-(et.st + epsilon)) - 1)) * (et.st/(et.st + 
    epsilon)^2)) + (p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^(-1/(et.st + epsilon)) * (log((p1^(-(et.st + 
    epsilon)) + p2^(-(et.st + epsilon)) - 1)) * (et.st/(et.st + 
    epsilon)^2 - et.st * (2 * (et.st * (et.st + 
    epsilon)))/((et.st + epsilon)^2)^2) - (p2^(-(et.st + 
    epsilon)) * (lo.p2 * et.st) + p1^(-(et.st + 
    epsilon)) * (lo.p1 * et.st))/(p1^(-(et.st + 
    epsilon)) + p2^(-(et.st + epsilon)) - 1) * (et.st/(et.st + 
    epsilon)^2)) - (((p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^((-1/(et.st + epsilon)) - 1) * (log((p1^(-(et.st + 
    epsilon)) + p2^(-(et.st + epsilon)) - 1)) * (et.st/(et.st + 
    epsilon)^2)) - (p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^(((-1/(et.st + epsilon)) - 1) - 1) * 
    (((-1/(et.st + epsilon)) - 1) * (p2^(-(et.st + 
        epsilon)) * (lo.p2 * et.st) + p1^(-(et.st + 
        epsilon)) * (lo.p1 * et.st)))) * ((-1/(et.st + 
    epsilon)) * (p2^(-(et.st + epsilon)) * (lo.p2 * 
    et.st) + p1^(-(et.st + epsilon)) * (lo.p1 * 
    et.st))) + (p1^(-(et.st + epsilon)) + p2^(-(et.st + 
    epsilon)) - 1)^((-1/(et.st + epsilon)) - 1) * (et.st/(et.st + 
    epsilon)^2 * (p2^(-(et.st + epsilon)) * (lo.p2 * 
    et.st) + p1^(-(et.st + epsilon)) * (lo.p1 * 
    et.st)) + (-1/(et.st + epsilon)) * (p2^(-(et.st + 
    epsilon)) * (lo.p2 * et.st) - p2^(-(et.st + 
    epsilon)) * (lo.p2 * et.st) * (lo.p2 * et.st) + 
    (p1^(-(et.st + epsilon)) * (lo.p1 * et.st) - 
        p1^(-(et.st + epsilon)) * (lo.p1 * et.st) * 
            (lo.p1 * et.st)))))

bit1.th2ATE <- ((p1^(-teta) + p2^(-teta) - 1)^(-1/teta) * (log((p1^(-teta) + 
    p2^(-teta) - 1)) * (1/teta^2)) - (p1^(-teta) + p2^(-teta) - 
    1)^((-1/teta) - 1) * ((-1/teta) * (p2^(-teta) * log(p2) + 
    p1^(-teta) * log(p1)))) * (log((p1^(-teta) + p2^(-teta) - 
    1)) * (1/teta^2)) - (p1^(-teta) + p2^(-teta) - 1)^(-1/teta) * 
    (log((p1^(-teta) + p2^(-teta) - 1)) * (2 * teta/(teta^2)^2) + 
        (p2^(-teta) * log(p2) + p1^(-teta) * log(p1))/(p1^(-teta) + 
            p2^(-teta) - 1) * (1/teta^2)) - (((p1^(-teta) + p2^(-teta) - 
    1)^((-1/teta) - 1) * (log((p1^(-teta) + p2^(-teta) - 1)) * 
    (1/teta^2)) - (p1^(-teta) + p2^(-teta) - 1)^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (p2^(-teta) * log(p2) + p1^(-teta) * 
    log(p1)))) * ((-1/teta) * (p2^(-teta) * lo.p2 + p1^(-teta) * 
    lo.p1)) + (p1^(-teta) + p2^(-teta) - 1)^((-1/teta) - 1) * 
    (1/teta^2 * (p2^(-teta) * lo.p2 + p1^(-teta) * lo.p1) - (-1/teta) * 
        (p1^(-teta) * log(p1) * lo.p1 + p2^(-teta) * log(p2) * 
            lo.p2)))


}




 
if(BivD=="C90"){

lo.p2 <- log(p2) 
et.st <- exp(teta.st)

c.copula.be1 <- ((1 - p1)^(-(-teta)) + p2^(-(-teta)) - 1)^((-1/(-teta)) - 1) * 
    ((-1/(-teta)) * ((1 - p1)^((-(-teta)) - 1) * (-(-teta))))
  c.copula.be2 <- 1 - ((1 - p1)^(-(-teta)) + p2^(-(-teta)) - 1)^((-1/(-teta)) - 
    1) * ((-1/(-teta)) * (p2^((-(-teta)) - 1) * (-(-teta))))
  c.copula.theta <- (-(((1 - p1)^teta + p2^teta - 1)^((1/teta) - 1) * ((1/teta) * 
    ((1 - p1)^teta * log((1 - p1)) + p2^teta * lo.p2)) - ((1 - 
    p1)^teta + p2^teta - 1)^(1/teta) * (log(((1 - p1)^teta + 
    p2^teta - 1)) * (1/teta^2))))*(-et.st)

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
     1) * ((1 - p1)^teta * log((1 - p1)) + p2^teta * lo.p2)) - 
     ((1 - p1)^teta + p2^teta - 1)^((1/teta) - 1) * (log(((1 - 
         p1)^teta + p2^teta - 1)) * (1/teta^2))) * ((1/teta) * 
     ((1 - p1)^(teta - 1) * teta)) + ((1 - p1)^teta + p2^teta - 
     1)^((1/teta) - 1) * ((1/teta) * ((1 - p1)^(teta - 1) * log((1 - 
     p1)) * teta + (1 - p1)^(teta - 1)) - 1/teta^2 * ((1 - p1)^(teta - 
     1) * teta)))*(-et.st)
   
 c.copula2.be2th <- (-((((1 - p1)^teta + p2^teta - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
     1) * ((1 - p1)^teta * log((1 - p1)) + p2^teta * lo.p2)) - 
     ((1 - p1)^teta + p2^teta - 1)^((1/teta) - 1) * (log(((1 - 
         p1)^teta + p2^teta - 1)) * (1/teta^2))) * ((1/teta) * 
     (p2^(teta - 1) * teta)) + ((1 - p1)^teta + p2^teta - 1)^((1/teta) - 
     1) * ((1/teta) * (p2^(teta - 1) * lo.p2 * teta + p2^(teta - 
     1)) - 1/teta^2 * (p2^(teta - 1) * teta))))*(-et.st)
 
 bit1.th2 <- -((((1 - p1)^(-(et.st + epsilon)) + p2^(-(et.st + 
     epsilon)) - 1)^(1/(-(et.st + epsilon))) * (log(((1 - 
     p1)^(-(et.st + epsilon)) + p2^(-(et.st + epsilon)) - 
     1)) * (et.st/(-(et.st + epsilon))^2)) - ((1 - 
     p1)^(-(et.st + epsilon)) + p2^(-(et.st + epsilon)) - 
     1)^((1/(-(et.st + epsilon))) - 1) * ((1/(-(et.st + 
     epsilon))) * (p2^(-(et.st + epsilon)) * (lo.p2 * 
     et.st) + (1 - p1)^(-(et.st + epsilon)) * (log((1 - 
     p1)) * et.st)))) * (log(((1 - p1)^(-(et.st + 
     epsilon)) + p2^(-(et.st + epsilon)) - 1)) * (et.st/(-(et.st + 
     epsilon))^2)) + ((1 - p1)^(-(et.st + epsilon)) + p2^(-(et.st + 
     epsilon)) - 1)^(1/(-(et.st + epsilon))) * (log(((1 - 
     p1)^(-(et.st + epsilon)) + p2^(-(et.st + epsilon)) - 
     1)) * (et.st/(-(et.st + epsilon))^2 + et.st * 
     (2 * (et.st * (-(et.st + epsilon))))/((-(et.st + 
     epsilon))^2)^2) - (p2^(-(et.st + epsilon)) * (lo.p2 * 
     et.st) + (1 - p1)^(-(et.st + epsilon)) * (log((1 - 
     p1)) * et.st))/((1 - p1)^(-(et.st + epsilon)) + 
     p2^(-(et.st + epsilon)) - 1) * (et.st/(-(et.st + 
     epsilon))^2)) - ((((1 - p1)^(-(et.st + epsilon)) + 
     p2^(-(et.st + epsilon)) - 1)^((1/(-(et.st + 
     epsilon))) - 1) * (log(((1 - p1)^(-(et.st + epsilon)) + 
     p2^(-(et.st + epsilon)) - 1)) * (et.st/(-(et.st + 
     epsilon))^2)) - ((1 - p1)^(-(et.st + epsilon)) + p2^(-(et.st + 
     epsilon)) - 1)^(((1/(-(et.st + epsilon))) - 1) - 1) * 
     (((1/(-(et.st + epsilon))) - 1) * (p2^(-(et.st + 
         epsilon)) * (lo.p2 * et.st) + (1 - p1)^(-(et.st + 
         epsilon)) * (log((1 - p1)) * et.st)))) * ((1/(-(et.st + 
     epsilon))) * (p2^(-(et.st + epsilon)) * (lo.p2 * 
     et.st) + (1 - p1)^(-(et.st + epsilon)) * (log((1 - 
     p1)) * et.st))) + ((1 - p1)^(-(et.st + epsilon)) + 
     p2^(-(et.st + epsilon)) - 1)^((1/(-(et.st + 
     epsilon))) - 1) * (et.st/(-(et.st + epsilon))^2 * 
     (p2^(-(et.st + epsilon)) * (lo.p2 * et.st) + 
         (1 - p1)^(-(et.st + epsilon)) * (log((1 - p1)) * 
             et.st)) + (1/(-(et.st + epsilon))) * 
     (p2^(-(et.st + epsilon)) * (lo.p2 * et.st) - 
         p2^(-(et.st + epsilon)) * (lo.p2 * et.st) * 
             (lo.p2 * et.st) + ((1 - p1)^(-(et.st + 
         epsilon)) * (log((1 - p1)) * et.st) - (1 - p1)^(-(et.st + 
         epsilon)) * (log((1 - p1)) * et.st) * (log((1 - 
         p1)) * et.st))))))
 

 bit1.th2ATE <- -((((1 - p1)^teta + p2^teta - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
    1) * ((1 - p1)^teta * log((1 - p1)) + p2^teta * log(p2))) - 
    ((1 - p1)^teta + p2^teta - 1)^((1/teta) - 1) * (log(((1 - 
        p1)^teta + p2^teta - 1)) * (1/teta^2))) * ((1/teta) * 
    ((1 - p1)^teta * log((1 - p1)) + p2^teta * lo.p2)) + ((1 - 
    p1)^teta + p2^teta - 1)^((1/teta) - 1) * ((1/teta) * ((1 - 
    p1)^teta * log((1 - p1)) * log((1 - p1)) + p2^teta * log(p2) * 
    lo.p2) - 1/teta^2 * ((1 - p1)^teta * log((1 - p1)) + p2^teta * 
    lo.p2)) - ((((1 - p1)^teta + p2^teta - 1)^((1/teta) - 1) * 
    ((1/teta) * ((1 - p1)^teta * log((1 - p1)) + p2^teta * log(p2))) - 
    ((1 - p1)^teta + p2^teta - 1)^(1/teta) * (log(((1 - p1)^teta + 
        p2^teta - 1)) * (1/teta^2))) * (log(((1 - p1)^teta + 
    p2^teta - 1)) * (1/teta^2)) + ((1 - p1)^teta + p2^teta - 
    1)^(1/teta) * (((1 - p1)^teta * log((1 - p1)) + p2^teta * 
    log(p2))/((1 - p1)^teta + p2^teta - 1) * (1/teta^2) - log(((1 - 
    p1)^teta + p2^teta - 1)) * (2 * teta/(teta^2)^2)))) 
 


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


bit1.th2ATE <- (((1 - p1)^(-teta) + (1 - p2)^(-teta) - 1)^(-1/teta) * (log(((1 - 
    p1)^(-teta) + (1 - p2)^(-teta) - 1)) * (1/teta^2)) - ((1 - 
    p1)^(-teta) + (1 - p2)^(-teta) - 1)^((-1/teta) - 1) * ((-1/teta) * 
    ((1 - p2)^(-teta) * log((1 - p2)) + (1 - p1)^(-teta) * log((1 - 
        p1))))) * (log(((1 - p1)^(-teta) + (1 - p2)^(-teta) - 
    1)) * (1/teta^2)) - ((1 - p1)^(-teta) + (1 - p2)^(-teta) - 
    1)^(-1/teta) * (log(((1 - p1)^(-teta) + (1 - p2)^(-teta) - 
    1)) * (2 * teta/(teta^2)^2) + ((1 - p2)^(-teta) * log((1 - 
    p2)) + (1 - p1)^(-teta) * log((1 - p1)))/((1 - p1)^(-teta) + 
    (1 - p2)^(-teta) - 1) * (1/teta^2)) - ((((1 - p1)^(-teta) + 
    (1 - p2)^(-teta) - 1)^((-1/teta) - 1) * (log(((1 - p1)^(-teta) + 
    (1 - p2)^(-teta) - 1)) * (1/teta^2)) - ((1 - p1)^(-teta) + 
    (1 - p2)^(-teta) - 1)^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * ((1 - p2)^(-teta) * log((1 - p2)) + (1 - p1)^(-teta) * 
    log((1 - p1))))) * ((-1/teta) * ((1 - p2)^(-teta) * log((1 - 
    p2)) + (1 - p1)^(-teta) * log((1 - p1)))) + ((1 - p1)^(-teta) + 
    (1 - p2)^(-teta) - 1)^((-1/teta) - 1) * (1/teta^2 * ((1 - 
    p2)^(-teta) * log((1 - p2)) + (1 - p1)^(-teta) * log((1 - 
    p1))) - (-1/teta) * ((1 - p1)^(-teta) * log((1 - p1)) * log((1 - 
    p1)) + (1 - p2)^(-teta) * log((1 - p2)) * log((1 - p2)))))
 
 


}




if(BivD=="C270"){

lo.p1 <- log(p1)

c.copula.be1 <- 1 - (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * ((1/teta) * 
    (p1^((teta) - 1) * (teta)))

  c.copula.be2 <- (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * ((1/teta) * 
    ((1 - p2)^((teta) - 1) * (teta)))




  c.copula.theta <- (-((p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * ((1/teta) * 
    (p1^(teta) * lo.p1 + (1 - p2)^(teta) * log((1 - p2)))) - 
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
     1) * (p1^(teta) * lo.p1 + (1 - p2)^(teta) * log((1 - p2)))) - 
     (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * (log((p1^(teta) + 
         (1 - p2)^(teta) - 1)) * (1/teta^2))) * ((1/teta) * (p1^((teta) - 
     1) * (teta))) + (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 
     1) * ((1/teta) * (p1^((teta) - 1) * lo.p1 * (teta) + p1^((teta) - 
     1)) - 1/teta^2 * (p1^((teta) - 1) * (teta)))))*(-(exp(teta.st)))
 
 
 c.copula2.be2th <- (((p1^(teta) + (1 - p2)^(teta) - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
     1) * (p1^(teta) * lo.p1 + (1 - p2)^(teta) * log((1 - p2)))) - 
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
     p1^(-(exp(teta.st) + epsilon)) * (lo.p1 * exp(teta.st))))) * 
     (log((p1^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
         epsilon)) - 1)) * (exp(teta.st)/(-(exp(teta.st) + epsilon))^2)) + 
     (p1^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
         epsilon)) - 1)^(1/(-(exp(teta.st) + epsilon))) * (log((p1^(-(exp(teta.st) + 
         epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)) * 
         (exp(teta.st)/(-(exp(teta.st) + epsilon))^2 + exp(teta.st) * 
             (2 * (exp(teta.st) * (-(exp(teta.st) + epsilon))))/((-(exp(teta.st) + 
             epsilon))^2)^2) - ((1 - p2)^(-(exp(teta.st) + epsilon)) * 
         (log((1 - p2)) * exp(teta.st)) + p1^(-(exp(teta.st) + 
         epsilon)) * (lo.p1 * exp(teta.st)))/(p1^(-(exp(teta.st) + 
         epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1) * 
         (exp(teta.st)/(-(exp(teta.st) + epsilon))^2)) - (((p1^(-(exp(teta.st) + 
     epsilon)) + (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)^((1/(-(exp(teta.st) + 
     epsilon))) - 1) * (log((p1^(-(exp(teta.st) + epsilon)) + 
     (1 - p2)^(-(exp(teta.st) + epsilon)) - 1)) * (exp(teta.st)/(-(exp(teta.st) + 
     epsilon))^2)) - (p1^(-(exp(teta.st) + epsilon)) + (1 - p2)^(-(exp(teta.st) + 
     epsilon)) - 1)^(((1/(-(exp(teta.st) + epsilon))) - 1) - 1) * 
     (((1/(-(exp(teta.st) + epsilon))) - 1) * ((1 - p2)^(-(exp(teta.st) + 
         epsilon)) * (log((1 - p2)) * exp(teta.st)) + p1^(-(exp(teta.st) + 
         epsilon)) * (lo.p1 * exp(teta.st))))) * ((1/(-(exp(teta.st) + 
     epsilon))) * ((1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - 
     p2)) * exp(teta.st)) + p1^(-(exp(teta.st) + epsilon)) * (lo.p1 * 
     exp(teta.st)))) + (p1^(-(exp(teta.st) + epsilon)) + (1 - 
     p2)^(-(exp(teta.st) + epsilon)) - 1)^((1/(-(exp(teta.st) + 
     epsilon))) - 1) * (exp(teta.st)/(-(exp(teta.st) + epsilon))^2 * 
     ((1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * 
         exp(teta.st)) + p1^(-(exp(teta.st) + epsilon)) * (lo.p1 * 
         exp(teta.st))) + (1/(-(exp(teta.st) + epsilon))) * ((1 - 
     p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * exp(teta.st)) - 
     (1 - p2)^(-(exp(teta.st) + epsilon)) * (log((1 - p2)) * exp(teta.st)) * 
         (log((1 - p2)) * exp(teta.st)) + (p1^(-(exp(teta.st) + 
     epsilon)) * (lo.p1 * exp(teta.st)) - p1^(-(exp(teta.st) + 
     epsilon)) * (lo.p1 * exp(teta.st)) * (lo.p1 * exp(teta.st)))))))
 
 bit1.th2ATE <- -(((p1^(teta) + (1 - p2)^(teta) - 1)^(((1/teta) - 1) - 1) * (((1/teta) - 
    1) * (p1^(teta) * log(p1) + (1 - p2)^(teta) * log((1 - p2)))) - 
    (p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 1) * (log((p1^(teta) + 
        (1 - p2)^(teta) - 1)) * (1/teta^2))) * ((1/teta) * (p1^(teta) * 
    lo.p1 + (1 - p2)^(teta) * log((1 - p2)))) + (p1^(teta) + 
    (1 - p2)^(teta) - 1)^((1/teta) - 1) * ((1/teta) * (p1^(teta) * 
    log(p1) * lo.p1 + (1 - p2)^(teta) * log((1 - p2)) * log((1 - 
    p2))) - 1/teta^2 * (p1^(teta) * lo.p1 + (1 - p2)^(teta) * 
    log((1 - p2)))) - (((p1^(teta) + (1 - p2)^(teta) - 1)^((1/teta) - 
    1) * ((1/teta) * (p1^(teta) * log(p1) + (1 - p2)^(teta) * 
    log((1 - p2)))) - (p1^(teta) + (1 - p2)^(teta) - 1)^(1/teta) * 
    (log((p1^(teta) + (1 - p2)^(teta) - 1)) * (1/teta^2))) * 
    (log((p1^(teta) + (1 - p2)^(teta) - 1)) * (1/teta^2)) + (p1^(teta) + 
    (1 - p2)^(teta) - 1)^(1/teta) * ((p1^(teta) * log(p1) + (1 - 
    p2)^(teta) * log((1 - p2)))/(p1^(teta) + (1 - p2)^(teta) - 
    1) * (1/teta^2) - log((p1^(teta) + (1 - p2)^(teta) - 1)) * 
    (2 * teta/(teta^2)^2))))
 
}














if(BivD=="F"){

# 1 - exp(-teta) = -expm1(-teta)
# recall log1p as well if needed. 

epcl <- -expm1(-teta) 


  c.copula.be1 <- -(-1/teta * (1/(epcl) * (exp(-teta * p1) * teta * (1 - 
    exp(-teta * p2)))/(1/(epcl) * ((epcl) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))))


 
  c.copula.be2 <- -(-1/teta * (1/(epcl) * ((1 - exp(-teta * p1)) * (exp(-teta * 
    p2) * teta))/(1/(epcl) * ((epcl) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))))



  c.copula.theta <- 1/teta^2 * log(1/(epcl) * ((epcl) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2)))) + -1/teta * ((1/(1 - 
    exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) - exp(-teta)/(epcl)^2 * ((epcl) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(epcl) * 
    ((epcl) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2)))))
        
c.copula2.be1 <-  -1/teta * (1/(epcl) * (exp(-teta * p1) * teta * teta * 
    (1 - exp(-teta * p2)))/(1/(epcl) * ((epcl) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) - 1/(1 - 
    exp(-teta)) * (exp(-teta * p1) * teta * (1 - exp(-teta * 
    p2))) * (1/(epcl) * (exp(-teta * p1) * teta * (1 - 
    exp(-teta * p2))))/(1/(epcl) * ((epcl) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))^2)
                  
 c.copula2.be2 <- -1/teta * (1/(epcl) * ((1 - exp(-teta * p1)) * (exp(-teta * 
    p2) * teta * teta))/(1/(epcl) * ((epcl) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) - 1/(1 - 
    exp(-teta)) * ((1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    teta)) * (1/(epcl) * ((1 - exp(-teta * p1)) * (exp(-teta * 
    p2) * teta)))/(1/(epcl) * ((epcl) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))^2)


c.copula2.be1be2 <- -(-1/teta * (1/(epcl) * (exp(-teta * p1) * teta * (exp(-teta * 
    p2) * teta))/(1/(epcl) * ((epcl) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2)))) + 1/(epcl) * 
    (exp(-teta * p1) * teta * (1 - exp(-teta * p2))) * (1/(1 - 
    exp(-teta)) * ((1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    teta)))/(1/(epcl) * ((epcl) - (1 - exp(-teta * 
    p1)) * (1 - exp(-teta * p2))))^2))

c.copula2.be1th <- -(1/teta^2 * (1/(epcl) * (exp(-teta * p1) * teta * 
    (1 - exp(-teta * p2)))/(1/(epcl) * ((epcl) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))) + -1/teta * 
    ((1/(epcl) * ((exp(-teta * p1) - exp(-teta * p1) * 
        p1 * teta) * (1 - exp(-teta * p2)) + exp(-teta * p1) * 
        teta * (exp(-teta * p2) * p2)) - exp(-teta)/(epcl)^2 * 
        (exp(-teta * p1) * teta * (1 - exp(-teta * p2))))/(1/(1 - 
        exp(-teta)) * ((epcl) - (1 - exp(-teta * p1)) * 
        (1 - exp(-teta * p2)))) - 1/(epcl) * (exp(-teta * 
        p1) * teta * (1 - exp(-teta * p2))) * (1/(epcl) * 
        (exp(-teta) - (exp(-teta * p1) * p1 * (1 - exp(-teta * 
            p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
            p2))) - exp(-teta)/(epcl)^2 * ((epcl) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - 
        exp(-teta)) * ((epcl) - (1 - exp(-teta * p1)) * 
        (1 - exp(-teta * p2))))^2))

c.copula2.be2th <- -(1/teta^2 * (1/(epcl) * ((1 - exp(-teta * p1)) * (exp(-teta * 
    p2) * teta))/(1/(epcl) * ((epcl) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))) + -1/teta * ((1/(1 - 
    exp(-teta)) * (exp(-teta * p1) * p1 * (exp(-teta * p2) * 
    teta) + (1 - exp(-teta * p1)) * (exp(-teta * p2) - exp(-teta * 
    p2) * p2 * teta)) - exp(-teta)/(epcl)^2 * ((1 - 
    exp(-teta * p1)) * (exp(-teta * p2) * teta)))/(1/(epcl) * 
    ((epcl) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2)))) - 1/(epcl) * ((1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * teta)) * (1/(epcl) * (exp(-teta) - 
    (exp(-teta * p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * 
        p1)) * (exp(-teta * p2) * p2))) - exp(-teta)/(epcl)^2 * 
    ((epcl) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))/(1/(epcl) * ((epcl) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))^2))


bit1.th2 <- bit1.th2ATE <- 1/teta^2 * ((1/(epcl) * (exp(-teta) - (exp(-teta * 
    p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * p2))) - exp(-teta)/(epcl)^2 * 
    ((epcl) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))/(1/(epcl) * ((epcl) - (1 - 
    exp(-teta * p1)) * (1 - exp(-teta * p2))))) - 2 * teta/(teta^2)^2 * 
    log(1/(epcl) * ((epcl) - (1 - exp(-teta * 
        p1)) * (1 - exp(-teta * p2)))) + (1/teta^2 * ((1/(1 - 
    exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) - exp(-teta)/(epcl)^2 * ((epcl) - 
    (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(epcl) * 
    ((epcl) - (1 - exp(-teta * p1)) * (1 - exp(-teta * 
        p2))))) - -1/teta * ((1/(epcl) * (exp(-teta) + 
    (exp(-teta * p1) * p1 * (exp(-teta * p2) * p2) - exp(-teta * 
        p1) * p1 * p1 * (1 - exp(-teta * p2)) + (exp(-teta * 
        p1) * p1 * (exp(-teta * p2) * p2) - (1 - exp(-teta * 
        p1)) * (exp(-teta * p2) * p2 * p2)))) + exp(-teta)/(1 - 
    exp(-teta))^2 * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
    exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * p2) * 
    p2))) + (exp(-teta)/(epcl)^2 * (exp(-teta) - (exp(-teta * 
    p1) * p1 * (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * 
    (exp(-teta * p2) * p2))) - (exp(-teta)/(epcl)^2 + 
    exp(-teta) * (2 * (exp(-teta) * (epcl)))/((1 - 
        exp(-teta))^2)^2) * ((epcl) - (1 - exp(-teta * 
    p1)) * (1 - exp(-teta * p2)))))/(1/(epcl) * ((1 - 
    exp(-teta)) - (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) + 
    (1/(epcl) * (exp(-teta) - (exp(-teta * p1) * p1 * 
        (1 - exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * 
        p2) * p2))) - exp(-teta)/(epcl)^2 * ((epcl) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2)))) * (1/(1 - 
        exp(-teta)) * (exp(-teta) - (exp(-teta * p1) * p1 * (1 - 
        exp(-teta * p2)) + (1 - exp(-teta * p1)) * (exp(-teta * 
        p2) * p2))) - exp(-teta)/(epcl)^2 * ((epcl) - 
        (1 - exp(-teta * p1)) * (1 - exp(-teta * p2))))/(1/(1 - 
        exp(-teta)) * ((epcl) - (1 - exp(-teta * p1)) * 
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

bit1.th2ATE <- -(exp(-((-log(p1))^teta + (-log(p2))^teta)^(1/teta)) * ((((-log(p1))^teta + 
    (-log(p2))^teta)^(((1/teta) - 1) - 1) * (((1/teta) - 1) * 
    ((-log(p1))^teta * log((-log(p1))) + (-log(p2))^teta * log((-log(p2))))) - 
    ((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 1) * (log(((-log(p1))^teta + 
        (-log(p2))^teta)) * (1/teta^2))) * ((1/teta) * ((-log(p1))^teta * 
    log((-log(p1))) + (-log(p2))^teta * log((-log(p2))))) + ((-log(p1))^teta + 
    (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p1))^teta * 
    log((-log(p1))) * log((-log(p1))) + (-log(p2))^teta * log((-log(p2))) * 
    log((-log(p2)))) - 1/teta^2 * ((-log(p1))^teta * log((-log(p1))) + 
    (-log(p2))^teta * log((-log(p2))))) - ((((-log(p1))^teta + 
    (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(p1))^teta * 
    log((-log(p1))) + (-log(p2))^teta * log((-log(p2))))) - ((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta) * (log(((-log(p1))^teta + (-log(p2))^teta)) * 
    (1/teta^2))) * (log(((-log(p1))^teta + (-log(p2))^teta)) * 
    (1/teta^2)) + ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
    (((-log(p1))^teta * log((-log(p1))) + (-log(p2))^teta * log((-log(p2))))/((-log(p1))^teta + 
        (-log(p2))^teta) * (1/teta^2) - log(((-log(p1))^teta + 
        (-log(p2))^teta)) * (2 * teta/(teta^2)^2)))) - exp(-((-log(p1))^teta + 
    (-log(p2))^teta)^(1/teta)) * (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 
    1) * ((1/teta) * ((-log(p1))^teta * log((-log(p1))) + (-log(p2))^teta * 
    log((-log(p2))))) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
    (log(((-log(p1))^teta + (-log(p2))^teta)) * (1/teta^2))) * 
    (((-log(p1))^teta + (-log(p2))^teta)^((1/teta) - 1) * ((1/teta) * 
        ((-log(p1))^teta * log((-log(p1))) + (-log(p2))^teta * 
            log((-log(p2))))) - ((-log(p1))^teta + (-log(p2))^teta)^(1/teta) * 
        (log(((-log(p1))^teta + (-log(p2))^teta)) * (1/teta^2))))

 
 

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

bit1.th2ATE <- exp(-((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta)) * 
    ((((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta) * 
        (log(((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))) * 
            (1/teta^2)) - ((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p2))^(-teta) * log((-log(p2))) + 
        (-log(1 - p1))^(-teta) * log((-log(1 - p1)))))) * (log(((-log(1 - 
        p1))^(-teta) + (-log(p2))^(-teta))) * (1/teta^2)) - ((-log(1 - 
        p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta) * (log(((-log(1 - 
        p1))^(-teta) + (-log(p2))^(-teta))) * (2 * teta/(teta^2)^2) + 
        ((-log(p2))^(-teta) * log((-log(p2))) + (-log(1 - p1))^(-teta) * 
            log((-log(1 - p1))))/((-log(1 - p1))^(-teta) + (-log(p2))^(-teta)) * 
            (1/teta^2)) - ((((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * (log(((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))) * 
        (1/teta^2)) - ((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(p2))^(-teta) * log((-log(p2))) + 
        (-log(1 - p1))^(-teta) * log((-log(1 - p1)))))) * ((-1/teta) * 
        ((-log(p2))^(-teta) * log((-log(p2))) + (-log(1 - p1))^(-teta) * 
            log((-log(1 - p1))))) + ((-log(1 - p1))^(-teta) + 
        (-log(p2))^(-teta))^((-1/teta) - 1) * (1/teta^2 * ((-log(p2))^(-teta) * 
        log((-log(p2))) + (-log(1 - p1))^(-teta) * log((-log(1 - 
        p1)))) - (-1/teta) * ((-log(1 - p1))^(-teta) * log((-log(1 - 
        p1))) * log((-log(1 - p1))) + (-log(p2))^(-teta) * log((-log(p2))) * 
        log((-log(p2))))))) - exp(-((-log(1 - p1))^(-teta) + 
    (-log(p2))^(-teta))^(-1/teta)) * (((-log(1 - p1))^(-teta) + 
    (-log(p2))^(-teta))^(-1/teta) * (log(((-log(1 - p1))^(-teta) + 
    (-log(p2))^(-teta))) * (1/teta^2)) - ((-log(1 - p1))^(-teta) + 
    (-log(p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((-log(p2))^(-teta) * 
    log((-log(p2))) + (-log(1 - p1))^(-teta) * log((-log(1 - 
    p1)))))) * (((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^(-1/teta) * 
    (log(((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))) * (1/teta^2)) - 
    ((-log(1 - p1))^(-teta) + (-log(p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(p2))^(-teta) * log((-log(p2))) + 
        (-log(1 - p1))^(-teta) * log((-log(1 - p1))))))

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

bit1.th2ATE <- -(exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
    ((((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * ((-log(1 - p1))^teta * log((-log(1 - 
        p1))) + (-log(1 - p2))^teta * log((-log(1 - p2))))) - 
        ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
            1) * (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * 
            (1/teta^2))) * ((1/teta) * ((-log(1 - p1))^teta * 
        log((-log(1 - p1))) + (-log(1 - p2))^teta * log((-log(1 - 
        p2))))) + ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
        1) * ((1/teta) * ((-log(1 - p1))^teta * log((-log(1 - 
        p1))) * log((-log(1 - p1))) + (-log(1 - p2))^teta * log((-log(1 - 
        p2))) * log((-log(1 - p2)))) - 1/teta^2 * ((-log(1 - 
        p1))^teta * log((-log(1 - p1))) + (-log(1 - p2))^teta * 
        log((-log(1 - p2))))) - ((((-log(1 - p1))^teta + (-log(1 - 
        p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - p1))^teta * 
        log((-log(1 - p1))) + (-log(1 - p2))^teta * log((-log(1 - 
        p2))))) - ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta) * 
        (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * (1/teta^2))) * 
        (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * (1/teta^2)) + 
        ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta) * 
            (((-log(1 - p1))^teta * log((-log(1 - p1))) + (-log(1 - 
                p2))^teta * log((-log(1 - p2))))/((-log(1 - p1))^teta + 
                (-log(1 - p2))^teta) * (1/teta^2) - log(((-log(1 - 
                p1))^teta + (-log(1 - p2))^teta)) * (2 * teta/(teta^2)^2)))) - 
    exp(-((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta)) * 
        (((-log(1 - p1))^teta + (-log(1 - p2))^teta)^((1/teta) - 
            1) * ((1/teta) * ((-log(1 - p1))^teta * log((-log(1 - 
            p1))) + (-log(1 - p2))^teta * log((-log(1 - p2))))) - 
            ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta) * 
                (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * 
                  (1/teta^2))) * (((-log(1 - p1))^teta + (-log(1 - 
        p2))^teta)^((1/teta) - 1) * ((1/teta) * ((-log(1 - p1))^teta * 
        log((-log(1 - p1))) + (-log(1 - p2))^teta * log((-log(1 - 
        p2))))) - ((-log(1 - p1))^teta + (-log(1 - p2))^teta)^(1/teta) * 
        (log(((-log(1 - p1))^teta + (-log(1 - p2))^teta)) * (1/teta^2))))
 

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


bit1.th2ATE <- exp(-((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta)) * 
    ((((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta) * 
        (log(((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))) * 
            (1/teta^2)) - ((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(1 - p2))^(-teta) * log((-log(1 - 
        p2))) + (-log(p1))^(-teta) * log((-log(p1)))))) * (log(((-log(p1))^(-teta) + 
        (-log(1 - p2))^(-teta))) * (1/teta^2)) - ((-log(p1))^(-teta) + 
        (-log(1 - p2))^(-teta))^(-1/teta) * (log(((-log(p1))^(-teta) + 
        (-log(1 - p2))^(-teta))) * (2 * teta/(teta^2)^2) + ((-log(1 - 
        p2))^(-teta) * log((-log(1 - p2))) + (-log(p1))^(-teta) * 
        log((-log(p1))))/((-log(p1))^(-teta) + (-log(1 - p2))^(-teta)) * 
        (1/teta^2)) - ((((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
        1) * (log(((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))) * 
        (1/teta^2)) - ((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((-log(1 - p2))^(-teta) * 
        log((-log(1 - p2))) + (-log(p1))^(-teta) * log((-log(p1)))))) * 
        ((-1/teta) * ((-log(1 - p2))^(-teta) * log((-log(1 - 
            p2))) + (-log(p1))^(-teta) * log((-log(p1))))) + 
        ((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
            1) * (1/teta^2 * ((-log(1 - p2))^(-teta) * log((-log(1 - 
            p2))) + (-log(p1))^(-teta) * log((-log(p1)))) - (-1/teta) * 
            ((-log(p1))^(-teta) * log((-log(p1))) * log((-log(p1))) + 
                (-log(1 - p2))^(-teta) * log((-log(1 - p2))) * 
                  log((-log(1 - p2))))))) - exp(-((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^(-1/teta)) * (((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^(-1/teta) * (log(((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))) * (1/teta^2)) - ((-log(p1))^(-teta) + 
    (-log(1 - p2))^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((-log(1 - 
    p2))^(-teta) * log((-log(1 - p2))) + (-log(p1))^(-teta) * 
    log((-log(p1)))))) * (((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^(-1/teta) * 
    (log(((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))) * (1/teta^2)) - 
    ((-log(p1))^(-teta) + (-log(1 - p2))^(-teta))^((-1/teta) - 
        1) * ((-1/teta) * ((-log(1 - p2))^(-teta) * log((-log(1 - 
        p2))) + (-log(p1))^(-teta) * log((-log(p1))))))
 

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
        
bit1.th2ATE <- -((((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - p1)^teta * log((1 - p1)) + 
    (1 - p2)^teta * log((1 - p2)) - ((1 - p1)^teta * log((1 - 
    p1)) * (1 - p2)^teta + (1 - p1)^teta * ((1 - p2)^teta * log((1 - 
    p2)))))) - ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^((1/teta) - 1) * (log(((1 - p1)^teta + (1 - 
    p2)^teta - (1 - p1)^teta * (1 - p2)^teta)) * (1/teta^2))) * 
    ((1/teta) * ((1 - p1)^teta * log((1 - p1)) + (1 - p2)^teta * 
        log((1 - p2)) - ((1 - p1)^teta * log((1 - p1)) * (1 - 
        p2)^teta + (1 - p1)^teta * ((1 - p2)^teta * log((1 - 
        p2)))))) + ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * 
    (1 - p2)^teta)^((1/teta) - 1) * ((1/teta) * ((1 - p1)^teta * 
    log((1 - p1)) * log((1 - p1)) + (1 - p2)^teta * log((1 - 
    p2)) * log((1 - p2)) - ((1 - p1)^teta * log((1 - p1)) * log((1 - 
    p1)) * (1 - p2)^teta + (1 - p1)^teta * log((1 - p1)) * ((1 - 
    p2)^teta * log((1 - p2))) + ((1 - p1)^teta * log((1 - p1)) * 
    ((1 - p2)^teta * log((1 - p2))) + (1 - p1)^teta * ((1 - p2)^teta * 
    log((1 - p2)) * log((1 - p2)))))) - 1/teta^2 * ((1 - p1)^teta * 
    log((1 - p1)) + (1 - p2)^teta * log((1 - p2)) - ((1 - p1)^teta * 
    log((1 - p1)) * (1 - p2)^teta + (1 - p1)^teta * ((1 - p2)^teta * 
    log((1 - p2)))))) - ((((1 - p1)^teta + (1 - p2)^teta - (1 - 
    p1)^teta * (1 - p2)^teta)^((1/teta) - 1) * ((1/teta) * ((1 - 
    p1)^teta * log((1 - p1)) + (1 - p2)^teta * log((1 - p2)) - 
    ((1 - p1)^teta * log((1 - p1)) * (1 - p2)^teta + (1 - p1)^teta * 
        ((1 - p2)^teta * log((1 - p2)))))) - ((1 - p1)^teta + 
    (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta) * 
    (log(((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - 
        p2)^teta)) * (1/teta^2))) * (log(((1 - p1)^teta + (1 - 
    p2)^teta - (1 - p1)^teta * (1 - p2)^teta)) * (1/teta^2)) + 
    ((1 - p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta)^(1/teta) * 
        (((1 - p1)^teta * log((1 - p1)) + (1 - p2)^teta * log((1 - 
            p2)) - ((1 - p1)^teta * log((1 - p1)) * (1 - p2)^teta + 
            (1 - p1)^teta * ((1 - p2)^teta * log((1 - p2)))))/((1 - 
            p1)^teta + (1 - p2)^teta - (1 - p1)^teta * (1 - p2)^teta) * 
            (1/teta^2) - log(((1 - p1)^teta + (1 - p2)^teta - 
            (1 - p1)^teta * (1 - p2)^teta)) * (2 * teta/(teta^2)^2))))
 
 

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
        
bit1.th2ATE <- ((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^(-1/teta) * 
    (log((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))) * 
        (1/teta^2)) - (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * 
    (1 - p2)^(-teta))^((-1/teta) - 1) * ((-1/teta) * ((1 - p2)^(-teta) * 
    log((1 - p2)) + p1^(-teta) * log(p1) - (p1^(-teta) * ((1 - 
    p2)^(-teta) * log((1 - p2))) + p1^(-teta) * log(p1) * (1 - 
    p2)^(-teta))))) * (log((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * 
    (1 - p2)^(-teta))) * (1/teta^2)) - (p1^(-teta) + (1 - p2)^(-teta) - 
    p1^(-teta) * (1 - p2)^(-teta))^(-1/teta) * (log((p1^(-teta) + 
    (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))) * (2 * 
    teta/(teta^2)^2) + ((1 - p2)^(-teta) * log((1 - p2)) + p1^(-teta) * 
    log(p1) - (p1^(-teta) * ((1 - p2)^(-teta) * log((1 - p2))) + 
    p1^(-teta) * log(p1) * (1 - p2)^(-teta)))/(p1^(-teta) + (1 - 
    p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta)) * (1/teta^2)) - 
    (((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
        1) * (log((p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * 
        (1 - p2)^(-teta))) * (1/teta^2)) - (p1^(-teta) + (1 - 
        p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * ((1 - p2)^(-teta) * log((1 - 
        p2)) + p1^(-teta) * log(p1) - (p1^(-teta) * ((1 - p2)^(-teta) * 
        log((1 - p2))) + p1^(-teta) * log(p1) * (1 - p2)^(-teta))))) * 
        ((-1/teta) * ((1 - p2)^(-teta) * log((1 - p2)) + p1^(-teta) * 
            log(p1) - (p1^(-teta) * ((1 - p2)^(-teta) * log((1 - 
            p2))) + p1^(-teta) * log(p1) * (1 - p2)^(-teta)))) + 
        (p1^(-teta) + (1 - p2)^(-teta) - p1^(-teta) * (1 - p2)^(-teta))^((-1/teta) - 
            1) * (1/teta^2 * ((1 - p2)^(-teta) * log((1 - p2)) + 
            p1^(-teta) * log(p1) - (p1^(-teta) * ((1 - p2)^(-teta) * 
            log((1 - p2))) + p1^(-teta) * log(p1) * (1 - p2)^(-teta))) - 
            (-1/teta) * (p1^(-teta) * log(p1) * log(p1) + (1 - 
                p2)^(-teta) * log((1 - p2)) * log((1 - p2)) - 
                (p1^(-teta) * log(p1) * ((1 - p2)^(-teta) * log((1 - 
                  p2))) + p1^(-teta) * log(p1) * log(p1) * (1 - 
                  p2)^(-teta) + (p1^(-teta) * ((1 - p2)^(-teta) * 
                  log((1 - p2)) * log((1 - p2))) + p1^(-teta) * 
                  log(p1) * ((1 - p2)^(-teta) * log((1 - p2))))))))
 
 
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

bit1.th2ATE <- -(((p1^teta + p2^teta - p1^teta * p2^teta)^(((1/teta) - 1) - 
    1) * (((1/teta) - 1) * (p1^teta * log(p1) + p2^teta * log(p2) - 
    (p1^teta * log(p1) * p2^teta + p1^teta * (p2^teta * log(p2))))) - 
    (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 1) * 
        (log((p1^teta + p2^teta - p1^teta * p2^teta)) * (1/teta^2))) * 
    ((1/teta) * (p1^teta * log(p1) + p2^teta * log(p2) - (p1^teta * 
        log(p1) * p2^teta + p1^teta * (p2^teta * log(p2))))) + 
    (p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 1) * 
        ((1/teta) * (p1^teta * log(p1) * log(p1) + p2^teta * 
            log(p2) * log(p2) - (p1^teta * log(p1) * log(p1) * 
            p2^teta + p1^teta * log(p1) * (p2^teta * log(p2)) + 
            (p1^teta * log(p1) * (p2^teta * log(p2)) + p1^teta * 
                (p2^teta * log(p2) * log(p2))))) - 1/teta^2 * 
            (p1^teta * log(p1) + p2^teta * log(p2) - (p1^teta * 
                log(p1) * p2^teta + p1^teta * (p2^teta * log(p2))))) - 
    (((p1^teta + p2^teta - p1^teta * p2^teta)^((1/teta) - 1) * 
        ((1/teta) * (p1^teta * log(p1) + p2^teta * log(p2) - 
            (p1^teta * log(p1) * p2^teta + p1^teta * (p2^teta * 
                log(p2))))) - (p1^teta + p2^teta - p1^teta * 
        p2^teta)^(1/teta) * (log((p1^teta + p2^teta - p1^teta * 
        p2^teta)) * (1/teta^2))) * (log((p1^teta + p2^teta - 
        p1^teta * p2^teta)) * (1/teta^2)) + (p1^teta + p2^teta - 
        p1^teta * p2^teta)^(1/teta) * ((p1^teta * log(p1) + p2^teta * 
        log(p2) - (p1^teta * log(p1) * p2^teta + p1^teta * (p2^teta * 
        log(p2))))/(p1^teta + p2^teta - p1^teta * p2^teta) * 
        (1/teta^2) - log((p1^teta + p2^teta - p1^teta * p2^teta)) * 
        (2 * teta/(teta^2)^2))))
 
 
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
        
bit1.th2ATE <- (((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^(-1/teta) * 
    (log(((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * 
        p2^(-teta))) * (1/teta^2)) - ((1 - p1)^(-teta) + p2^(-teta) - 
    (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 1) * ((-1/teta) * 
    (p2^(-teta) * log(p2) + (1 - p1)^(-teta) * log((1 - p1)) - 
        ((1 - p1)^(-teta) * (p2^(-teta) * log(p2)) + (1 - p1)^(-teta) * 
            log((1 - p1)) * p2^(-teta))))) * (log(((1 - p1)^(-teta) + 
    p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))) * (1/teta^2)) - 
    ((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^(-1/teta) * 
        (log(((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * 
            p2^(-teta))) * (2 * teta/(teta^2)^2) + (p2^(-teta) * 
            log(p2) + (1 - p1)^(-teta) * log((1 - p1)) - ((1 - 
            p1)^(-teta) * (p2^(-teta) * log(p2)) + (1 - p1)^(-teta) * 
            log((1 - p1)) * p2^(-teta)))/((1 - p1)^(-teta) + 
            p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta)) * (1/teta^2)) - 
    ((((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 
        1) * (log(((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * 
        p2^(-teta))) * (1/teta^2)) - ((1 - p1)^(-teta) + p2^(-teta) - 
        (1 - p1)^(-teta) * p2^(-teta))^(((-1/teta) - 1) - 1) * 
        (((-1/teta) - 1) * (p2^(-teta) * log(p2) + (1 - p1)^(-teta) * 
            log((1 - p1)) - ((1 - p1)^(-teta) * (p2^(-teta) * 
            log(p2)) + (1 - p1)^(-teta) * log((1 - p1)) * p2^(-teta))))) * 
        ((-1/teta) * (p2^(-teta) * log(p2) + (1 - p1)^(-teta) * 
            log((1 - p1)) - ((1 - p1)^(-teta) * (p2^(-teta) * 
            log(p2)) + (1 - p1)^(-teta) * log((1 - p1)) * p2^(-teta)))) + 
        ((1 - p1)^(-teta) + p2^(-teta) - (1 - p1)^(-teta) * p2^(-teta))^((-1/teta) - 
            1) * (1/teta^2 * (p2^(-teta) * log(p2) + (1 - p1)^(-teta) * 
            log((1 - p1)) - ((1 - p1)^(-teta) * (p2^(-teta) * 
            log(p2)) + (1 - p1)^(-teta) * log((1 - p1)) * p2^(-teta))) - 
            (-1/teta) * ((1 - p1)^(-teta) * log((1 - p1)) * log((1 - 
                p1)) + p2^(-teta) * log(p2) * log(p2) - ((1 - 
                p1)^(-teta) * log((1 - p1)) * (p2^(-teta) * log(p2)) + 
                (1 - p1)^(-teta) * log((1 - p1)) * log((1 - p1)) * 
                  p2^(-teta) + ((1 - p1)^(-teta) * (p2^(-teta) * 
                log(p2) * log(p2)) + (1 - p1)^(-teta) * log((1 - 
                p1)) * (p2^(-teta) * log(p2)))))))
 
 

}








if(PL=="PP"){

p1.1 <- pnorm(eta1)
d1.1 <- dnorm(eta1)
p2.2 <- pnorm(eta2)
d2.2 <- dnorm(eta2)

der.d.n1.be1 <- p1.1^((xi1 - 1) - 1) * ((xi1 - 1) * d1.1) * 
     (xi1 * d1.1) - p1.1^(xi1 - 1) * (xi1 * 
     (eta1 * d1.1)) 
     
der.d.n2.be2 <- p2.2^((xi2 - 1) - 1) * ((xi2 - 1) * d2.2) * 
     (xi2 * d2.2) - p2.2^(xi2 - 1) * (xi2 * 
     (eta2 * d2.2)) 


if(eqPL=="both"){

der.p1.lambda1 <- (p1 * log(p1^(1/xi1)))*exp(xi1.st)  

der.p2.lambda2 <- (p2 * log(p2^(1/xi2)))*exp(xi2.st)   
        
     
der.der.p1.lam1.der.p1 <- (log(p1^(1/xi1)) + p1 * (p1^((1/xi1) - 1) * (1/xi1)/p1^(1/xi1))) * 
    exp(xi1.st)
    
der.der.p2.lam2.der.p2 <-  (log(p2^(1/xi2)) + p2 * (p2^((1/xi2) - 1) * (1/xi2)/p2^(1/xi2))) * 
    exp(xi2.st)    
    
der2.p1.lambda1 <- p1.1^(exp(xi1.st) + epsilon) * (log(p1.1) * 
     exp(xi1.st)) * (log(p1.1) * exp(xi1.st)) + 
     p1.1^(exp(xi1.st) + epsilon) * (log(p1.1) * 
         exp(xi1.st))
         
der2.p2.lambda2 <- p2.2^(exp(xi2.st) + epsilon) * (log(p2.2) * 
 	    exp(xi2.st)) * (log(p2.2) * exp(xi2.st)) + 
 	    p2.2^(exp(xi2.st) + epsilon) * (log(p2.2) * 
         exp(xi2.st))
         
der.d.n1.lambda1 <- (p1.1^(xi1 - 1) * log(p1.1) * (xi1 * d1.1) + 
     p1.1^(xi1 - 1) * d1.1)*exp(xi1.st) 
     
der.d.n2.lambda2 <- (p2.2^(xi2 - 1) * log(p2.2) * (xi2 * d2.2) + 
     p2.2^(xi2 - 1) * d2.2)*exp(xi2.st)     
}

if(eqPL=="first"){
der.p1.lambda1 <- (p1 * log(p1^(1/xi1)))*exp(xi1.st)  

der.der.p1.lam1.der.p1 <- (log(p1^(1/xi1)) + p1 * (p1^((1/xi1) - 1) * (1/xi1)/p1^(1/xi1))) * 
    exp(xi1.st)
      
der2.p1.lambda1 <- p1.1^(exp(xi1.st) + epsilon) * (log(p1.1) * 
     exp(xi1.st)) * (log(p1.1) * exp(xi1.st)) + 
     p1.1^(exp(xi1.st) + epsilon) * (log(p1.1) * 
         exp(xi1.st))
           
der.d.n1.lambda1 <- (p1.1^(xi1 - 1) * log(p1.1) * (xi1 * d1.1) + 
     p1.1^(xi1 - 1) * d1.1)*exp(xi1.st) 
     
   
}

if(eqPL=="second"){

der.p2.lambda2 <- (p2 * log(p2^(1/xi2)))*exp(xi2.st)   
         
der.der.p2.lam2.der.p2 <-  (log(p2^(1/xi2)) + p2 * (p2^((1/xi2) - 1) * (1/xi2)/p2^(1/xi2))) * 
    exp(xi2.st)    
    
der2.p2.lambda2 <- p2.2^(exp(xi2.st) + epsilon) * (log(p2.2) * 
 	    exp(xi2.st)) * (log(p2.2) * exp(xi2.st)) + 
 	    p2.2^(exp(xi2.st) + epsilon) * (log(p2.2) * 
         exp(xi2.st))
             
der.d.n2.lambda2 <- (p2.2^(xi2 - 1) * log(p2.2) * (xi2 * d2.2) + 
     p2.2^(xi2 - 1) * d2.2)*exp(xi2.st)     
}

}






if(PL=="RPP"){

p1.1 <- pnorm(-eta1)
d1.1 <- dnorm(-eta1)
p2.2 <- pnorm(-eta2)
d2.2 <- dnorm(-eta2)

der.d.n1.be1 <- -(p1.1^(xi1 - 1) * (xi1 * (eta1 * d1.1)) + 
     p1.1^((xi1 - 1) - 1) * ((xi1 - 1) * d1.1) * 
         (xi1 * d1.1))
     
der.d.n2.be2 <- -(p2.2^(xi2 - 1) * (xi2 * (eta2 * d2.2)) + 
     p2.2^((xi2 - 1) - 1) * ((xi2 - 1) * d2.2) * 
         (xi2 * d2.2))  



if(eqPL=="both"){

der.p1.lambda1 <- (-((1-p1)*log((1-p1)^(1/xi1))))*exp(xi1.st)  
   
der.p2.lambda2 <- (-((1-p2)*log((1-p2)^(1/xi2))))*exp(xi2.st)  
                      
der.der.p1.lam1.der.p1 <- ((1 - p1) * ((1 - p1)^((1/xi1) - 1) * (1/xi1)/(1 - p1)^(1/xi1)) + 
    log((1 - p1)^(1/xi1))) * exp(xi1.st)#

der.der.p2.lam2.der.p2 <-  ((1 - p2) * ((1 - p2)^((1/xi2) - 1) * (1/xi2)/(1 - p2)^(1/xi2)) + 
    log((1 - p2)^(1/xi2))) * exp(xi2.st)

der2.p1.lambda1 <--(p1.1^(exp(xi1.st) + epsilon) * (log(p1.1) * 
    exp(xi1.st)) * (log(p1.1) * exp(xi1.st)) + 
    p1.1^(exp(xi1.st) + epsilon) * (log(p1.1) * 
        exp(xi1.st)))

der2.p2.lambda2 <- -(p2.2^(exp(xi2.st) + epsilon) * (log(p2.2) * 
    exp(xi2.st)) * (log(p2.2) * exp(xi2.st)) + 
    p2.2^(exp(xi2.st) + epsilon) * (log(p2.2) * 
        exp(xi2.st)))
   
der.d.n1.lambda1 <- (p1.1^(xi1 - 1) * log(p1.1) * (xi1 * d1.1) + 
    p1.1^(xi1 - 1) * d1.1)*exp(xi1.st)

der.d.n2.lambda2 <- (p2.2^(xi2 - 1) * log(p2.2) * (xi2 * d2.2) + 
    p2.2^(xi2 - 1) * d2.2)*exp(xi2.st)    
            
}

if(eqPL=="first"){

der.p1.lambda1 <- (-((1-p1)*log((1-p1)^(1/xi1))))*exp(xi1.st)  
                      
der.der.p1.lam1.der.p1 <- ((1 - p1) * ((1 - p1)^((1/xi1) - 1) * (1/xi1)/(1 - p1)^(1/xi1)) + 
    log((1 - p1)^(1/xi1))) * exp(xi1.st)

der2.p1.lambda1 <--(p1.1^(exp(xi1.st) + epsilon) * (log(p1.1) * 
    exp(xi1.st)) * (log(p1.1) * exp(xi1.st)) + 
    p1.1^(exp(xi1.st) + epsilon) * (log(p1.1) * 
        exp(xi1.st)))

der.d.n1.lambda1 <- (p1.1^(xi1 - 1) * log(p1.1) * (xi1 * d1.1) + 
    p1.1^(xi1 - 1) * d1.1)*exp(xi1.st)

   
            
}

if(eqPL=="second"){

der.p2.lambda2 <- (-((1-p2)*log((1-p2)^(1/xi2))))*exp(xi2.st)  
                      
der.der.p2.lam2.der.p2 <-  ((1 - p2) * ((1 - p2)^((1/xi2) - 1) * (1/xi2)/(1 - p2)^(1/xi2)) + 
    log((1 - p2)^(1/xi2))) * exp(xi2.st)

der2.p2.lambda2 <- -(p2.2^(exp(xi2.st) + epsilon) * (log(p2.2) * 
    exp(xi2.st)) * (log(p2.2) * exp(xi2.st)) + 
    p2.2^(exp(xi2.st) + epsilon) * (log(p2.2) * 
        exp(xi2.st)))
   
der.d.n2.lambda2 <- (p2.2^(xi2 - 1) * log(p2.2) * (xi2 * d2.2) + 
    p2.2^(xi2 - 1) * d2.2)*exp(xi2.st)    
            
}

}





if(PL=="SN"){


der.d.n1.be1 <- 2 * dnorm(eta1) * (dnorm(xi1 * eta1) * xi1) - 2 * (eta1 * 
    dnorm(eta1)) * pnorm(xi1 * eta1) 
   
der.d.n2.be2 <- 2 * dnorm(eta2) * (dnorm(xi2 * eta2) * xi2) - 2 * (eta2 * 
    dnorm(eta2)) * pnorm(xi2 * eta2) 


if(eqPL=="both"){

del1 <- -xi1/sqrt(1+xi1^2)

del2 <- -xi2/sqrt(1+xi2^2)

der.p1.lambda1 <- 2*dbinorm(eta1,0, cov12=del1)*(-(1/sqrt(1 + xi1^2) - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5))/sqrt(1 + xi1^2)^2))

der.p2.lambda2 <- 2*dbinorm(eta2,0, cov12=del2)*(-(1/sqrt(1 + xi2^2) - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5))/sqrt(1 + xi2^2)^2))

der.der.p1.lam1.der.p1 <- der.p1.lambda1*(-eta1/(1-del1^2))*((qsn(p1+eps,alpha=xi1)-qsn(p1,alpha=xi1))/eps)
der.der.p2.lam2.der.p2 <- der.p2.lambda2*(-eta2/(1-del2^2))*((qsn(p2+eps,alpha=xi2)-qsn(p2,alpha=xi2))/eps)

der2.p1.lambda1 <- 2 * ((1/(2 * pi * sqrt(1 - (-xi1/sqrt(1 + xi1^2))^2)) * (exp(-1/(2 * 
    (1 - (-xi1/sqrt(1 + xi1^2))^2)) * eta1^2) * (2 * (2 * ((1/sqrt(1 + 
    xi1^2) - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5))/sqrt(1 + 
    xi1^2)^2) * (-xi1/sqrt(1 + xi1^2))))/(2 * (1 - (-xi1/sqrt(1 + 
    xi1^2))^2))^2 * eta1^2)) - 2 * pi * (0.5 * (2 * ((1/sqrt(1 + 
    xi1^2) - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5))/sqrt(1 + 
    xi1^2)^2) * (-xi1/sqrt(1 + xi1^2))) * (1 - (-xi1/sqrt(1 + 
    xi1^2))^2)^-0.5))/(2 * pi * sqrt(1 - (-xi1/sqrt(1 + xi1^2))^2))^2 * 
    exp(-1/(2 * (1 - (-xi1/sqrt(1 + xi1^2))^2)) * eta1^2)) * 
    (-(1/sqrt(1 + xi1^2) - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5))/sqrt(1 + 
        xi1^2)^2)) + 1/(2 * pi * sqrt(1 - (-xi1/sqrt(1 + xi1^2))^2)) * 
    exp(-1/(2 * (1 - (-xi1/sqrt(1 + xi1^2))^2)) * eta1^2) * (0.5 * 
    (2 * xi1 * (1 + xi1^2)^-0.5)/sqrt(1 + xi1^2)^2 + (((0.5 * 
    (2 * xi1 * (1 + xi1^2)^-0.5)) + xi1 * (0.5 * (2 * (1 + xi1^2)^-0.5 - 
    2 * xi1 * ((1 + xi1^2)^-(0.5 + 1) * (0.5 * (2 * xi1))))))/sqrt(1 + 
    xi1^2)^2 - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5)) * (2 * 
    (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5) * sqrt(1 + xi1^2)))/(sqrt(1 + 
    xi1^2)^2)^2)))

  
der2.p2.lambda2 <- 2 * ((1/(2 * pi * sqrt(1 - (-xi2/sqrt(1 + xi2^2))^2)) * (exp(-1/(2 * 
    (1 - (-xi2/sqrt(1 + xi2^2))^2)) * eta2^2) * (2 * (2 * ((1/sqrt(1 + 
    xi2^2) - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5))/sqrt(1 + 
    xi2^2)^2) * (-xi2/sqrt(1 + xi2^2))))/(2 * (1 - (-xi2/sqrt(1 + 
    xi2^2))^2))^2 * eta2^2)) - 2 * pi * (0.5 * (2 * ((1/sqrt(1 + 
    xi2^2) - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5))/sqrt(1 + 
    xi2^2)^2) * (-xi2/sqrt(1 + xi2^2))) * (1 - (-xi2/sqrt(1 + 
    xi2^2))^2)^-0.5))/(2 * pi * sqrt(1 - (-xi2/sqrt(1 + xi2^2))^2))^2 * 
    exp(-1/(2 * (1 - (-xi2/sqrt(1 + xi2^2))^2)) * eta2^2)) * 
    (-(1/sqrt(1 + xi2^2) - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5))/sqrt(1 + 
        xi2^2)^2)) + 1/(2 * pi * sqrt(1 - (-xi2/sqrt(1 + xi2^2))^2)) * 
    exp(-1/(2 * (1 - (-xi2/sqrt(1 + xi2^2))^2)) * eta2^2) * (0.5 * 
    (2 * xi2 * (1 + xi2^2)^-0.5)/sqrt(1 + xi2^2)^2 + (((0.5 * 
    (2 * xi2 * (1 + xi2^2)^-0.5)) + xi2 * (0.5 * (2 * (1 + xi2^2)^-0.5 - 
    2 * xi2 * ((1 + xi2^2)^-(0.5 + 1) * (0.5 * (2 * xi2))))))/sqrt(1 + 
    xi2^2)^2 - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5)) * (2 * 
    (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5) * sqrt(1 + xi2^2)))/(sqrt(1 + 
    xi2^2)^2)^2)))

        

der.d.n1.lambda1 <- 2 * dnorm(eta1) * (dnorm(xi1 * eta1) * eta1) 
     
der.d.n2.lambda2 <- 2 * dnorm(eta2) * (dnorm(xi2 * eta2) * eta2)    


}




if(eqPL=="first"){


del1 <- -xi1/sqrt(1+xi1^2)


der.p1.lambda1 <- 2*dbinorm(eta1,0, cov12=del1)*(-(1/sqrt(1 + xi1^2) - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5))/sqrt(1 + xi1^2)^2))


der.der.p1.lam1.der.p1 <- der.p1.lambda1*(-eta1/(1-del1^2))*((qsn(p1+eps,alpha=xi1)-qsn(p1,alpha=xi1))/eps)


der2.p1.lambda1 <- 2 * ((1/(2 * pi * sqrt(1 - (-xi1/sqrt(1 + xi1^2))^2)) * (exp(-1/(2 * 
    (1 - (-xi1/sqrt(1 + xi1^2))^2)) * eta1^2) * (2 * (2 * ((1/sqrt(1 + 
    xi1^2) - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5))/sqrt(1 + 
    xi1^2)^2) * (-xi1/sqrt(1 + xi1^2))))/(2 * (1 - (-xi1/sqrt(1 + 
    xi1^2))^2))^2 * eta1^2)) - 2 * pi * (0.5 * (2 * ((1/sqrt(1 + 
    xi1^2) - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5))/sqrt(1 + 
    xi1^2)^2) * (-xi1/sqrt(1 + xi1^2))) * (1 - (-xi1/sqrt(1 + 
    xi1^2))^2)^-0.5))/(2 * pi * sqrt(1 - (-xi1/sqrt(1 + xi1^2))^2))^2 * 
    exp(-1/(2 * (1 - (-xi1/sqrt(1 + xi1^2))^2)) * eta1^2)) * 
    (-(1/sqrt(1 + xi1^2) - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5))/sqrt(1 + 
        xi1^2)^2)) + 1/(2 * pi * sqrt(1 - (-xi1/sqrt(1 + xi1^2))^2)) * 
    exp(-1/(2 * (1 - (-xi1/sqrt(1 + xi1^2))^2)) * eta1^2) * (0.5 * 
    (2 * xi1 * (1 + xi1^2)^-0.5)/sqrt(1 + xi1^2)^2 + (((0.5 * 
    (2 * xi1 * (1 + xi1^2)^-0.5)) + xi1 * (0.5 * (2 * (1 + xi1^2)^-0.5 - 
    2 * xi1 * ((1 + xi1^2)^-(0.5 + 1) * (0.5 * (2 * xi1))))))/sqrt(1 + 
    xi1^2)^2 - xi1 * (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5)) * (2 * 
    (0.5 * (2 * xi1 * (1 + xi1^2)^-0.5) * sqrt(1 + xi1^2)))/(sqrt(1 + 
    xi1^2)^2)^2)))

 
der.d.n1.lambda1 <- 2 * dnorm(eta1) * (dnorm(xi1 * eta1) * eta1) 
     
     
}




if(eqPL=="second"){



del2 <- -xi2/sqrt(1+xi2^2)


der.p2.lambda2 <- 2*dbinorm(eta2,0, cov12=del2)*(-(1/sqrt(1 + xi2^2) - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5))/sqrt(1 + xi2^2)^2))


der.der.p2.lam2.der.p2 <- der.p2.lambda2*(-eta2/(1-del2^2))*((qsn(p2+eps,alpha=xi2)-qsn(p2,alpha=xi2))/eps)

  
der2.p2.lambda2 <- 2 * ((1/(2 * pi * sqrt(1 - (-xi2/sqrt(1 + xi2^2))^2)) * (exp(-1/(2 * 
    (1 - (-xi2/sqrt(1 + xi2^2))^2)) * eta2^2) * (2 * (2 * ((1/sqrt(1 + 
    xi2^2) - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5))/sqrt(1 + 
    xi2^2)^2) * (-xi2/sqrt(1 + xi2^2))))/(2 * (1 - (-xi2/sqrt(1 + 
    xi2^2))^2))^2 * eta2^2)) - 2 * pi * (0.5 * (2 * ((1/sqrt(1 + 
    xi2^2) - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5))/sqrt(1 + 
    xi2^2)^2) * (-xi2/sqrt(1 + xi2^2))) * (1 - (-xi2/sqrt(1 + 
    xi2^2))^2)^-0.5))/(2 * pi * sqrt(1 - (-xi2/sqrt(1 + xi2^2))^2))^2 * 
    exp(-1/(2 * (1 - (-xi2/sqrt(1 + xi2^2))^2)) * eta2^2)) * 
    (-(1/sqrt(1 + xi2^2) - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5))/sqrt(1 + 
        xi2^2)^2)) + 1/(2 * pi * sqrt(1 - (-xi2/sqrt(1 + xi2^2))^2)) * 
    exp(-1/(2 * (1 - (-xi2/sqrt(1 + xi2^2))^2)) * eta2^2) * (0.5 * 
    (2 * xi2 * (1 + xi2^2)^-0.5)/sqrt(1 + xi2^2)^2 + (((0.5 * 
    (2 * xi2 * (1 + xi2^2)^-0.5)) + xi2 * (0.5 * (2 * (1 + xi2^2)^-0.5 - 
    2 * xi2 * ((1 + xi2^2)^-(0.5 + 1) * (0.5 * (2 * xi2))))))/sqrt(1 + 
    xi2^2)^2 - xi2 * (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5)) * (2 * 
    (0.5 * (2 * xi2 * (1 + xi2^2)^-0.5) * sqrt(1 + xi2^2)))/(sqrt(1 + 
    xi2^2)^2)^2)))

        
     
der.d.n2.lambda2 <- 2 * dnorm(eta2) * (dnorm(xi2 * eta2) * eta2)   



}



}





if(PL!="P"){

c.copula.lambda1 <- c.copula.be1*der.p1.lambda1 
c.copula.lambda2 <- c.copula.be2*der.p2.lambda2 


if(eqPL=="both"){

bit1.lambda1.2 <- c.copula2.be1*der.p1.lambda1^2+c.copula.be1*der2.p1.lambda1
bit1.lambda2.2 <- c.copula2.be2*der.p2.lambda2^2+c.copula.be2*der2.p2.lambda2       
c.copula2.be1lambda1 <- c.copula2.be1*der.p1.lambda1
c.copula2.be2lambda2 <- c.copula2.be2*der.p2.lambda2       
c.copula2.be1lambda2 <- c.copula2.be1be2*der.p2.lambda2       
c.copula2.be2lambda1 <- c.copula2.be1be2*der.p1.lambda1       
bit1.thlambda1 <- c.copula2.be1th* der.p1.lambda1
bit1.thlambda2 <- c.copula2.be2th* der.p2.lambda2
bit1.lambda1lambda2 <- c.copula2.be1be2*der.p1.lambda1*der.p2.lambda2

}

if(eqPL=="first"){
  
bit1.lambda1.2 <- c.copula2.be1*der.p1.lambda1^2+c.copula.be1*der2.p1.lambda1
c.copula2.be1lambda1 <- c.copula2.be1*der.p1.lambda1
c.copula2.be2lambda1 <- c.copula2.be1be2*der.p1.lambda1       
bit1.thlambda1 <- c.copula2.be1th* der.p1.lambda1


}

if(eqPL=="second"){

bit1.lambda2.2 <- c.copula2.be2*der.p2.lambda2^2+c.copula.be2*der2.p2.lambda2        
c.copula2.be2lambda2 <- c.copula2.be2*der.p2.lambda2        
c.copula2.be1lambda2 <- c.copula2.be1be2*der.p2.lambda2         
bit1.thlambda2 <- c.copula2.be2th* der.p2.lambda2

}

}



         list(c.copula.be1=c.copula.be1, c.copula.be2=c.copula.be2, c.copula.theta=c.copula.theta, 
              c.copula2.be1=c.copula2.be1, c.copula2.be2=c.copula2.be2, c.copula2.be1be2=c.copula2.be1be2,
              c.copula2.be1th=c.copula2.be1th,c.copula2.be2th=c.copula2.be2th,bit1.th2=bit1.th2,          
              c.copula.lambda1=c.copula.lambda1, c.copula.lambda2=c.copula.lambda2, bit1.lambda1.2=bit1.lambda1.2, 
              bit1.lambda2.2=bit1.lambda2.2, c.copula2.be1lambda1=c.copula2.be1lambda1, c.copula2.be2lambda2=c.copula2.be2lambda2,
              c.copula2.be1lambda2=c.copula2.be1lambda2, c.copula2.be2lambda1=c.copula2.be2lambda1, bit1.thlambda1=bit1.thlambda1, 
              bit1.thlambda2=bit1.thlambda2, bit1.lambda1lambda2=bit1.lambda1lambda2,
              der.p1.lambda1=der.p1.lambda1, der.p2.lambda2=der.p2.lambda2, der.d.n1.be1=der.d.n1.be1, der.d.n2.be2=der.d.n2.be2,
der.der.p1.lam1.der.p1=der.der.p1.lam1.der.p1, der.der.p2.lam2.der.p2=der.der.p2.lam2.der.p2, 
der2.p1.lambda1=der2.p1.lambda1, der2.p2.lambda2=der2.p2.lambda2, der.d.n1.lambda1=der.d.n1.lambda1,
der.d.n2.lambda2=der.d.n2.lambda2, bit1.th2ATE= bit1.th2ATE)     


}




     
























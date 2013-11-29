copgHsBB <- function(p1,p2,teta,teta.st,delta,delta.st,BivD){

epsilon <- 0 # .Machine$double.eps*10^6
 
if(BivD=="BB1.0"){

  c.copula.be1 <- -((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * ((p1^-teta - 1)^(delta - 1) * (delta * 
    (p1^-(teta + 1) * teta)))))))

 
  c.copula.be2 <- -((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * ((p2^-teta - 1)^(delta - 1) * (delta * 
    (p2^-(teta + 1) * teta)))))))

  c.copula.theta <- ((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^(-1/teta) * 
    (log((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))) * 
        (1/teta^2)) - (1 + ((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^(1/delta))^((-1/teta) - 1) * ((-1/teta) * (((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    ((p2^-teta - 1)^(delta - 1) * (delta * (p2^-teta * log(p2))) + 
        (p1^-teta - 1)^(delta - 1) * (delta * (p1^-teta * log(p1))))))))*exp(teta.st)


  c.copula.delta <- ((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * ((p1^-teta - 1)^delta * log((p1^-teta - 
    1)) + (p2^-teta - 1)^delta * log((p2^-teta - 1)))) - ((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)^(1/delta) * (log(((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)) * (1/delta^2)))))*exp(delta.st)


 c.copula2.be1 <-  (1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * ((p1^-teta - 1)^(delta - 1) * (delta * 
    (p1^-(teta + 1 + 1) * (teta + 1) * teta)) + (p1^-teta - 1)^((delta - 
    1) - 1) * ((delta - 1) * (p1^-(teta + 1) * teta)) * (delta * 
    (p1^-(teta + 1) * teta)))) + ((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((p1^-teta - 
    1)^(delta - 1) * (delta * (p1^-(teta + 1) * teta)))) * ((1/delta) * 
    ((p1^-teta - 1)^(delta - 1) * (delta * (p1^-(teta + 1) * 
        teta)))))) + (1 + ((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^(1/delta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * ((p1^-teta - 1)^(delta - 1) * (delta * 
    (p1^-(teta + 1) * teta)))))) * ((-1/teta) * (((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    ((p1^-teta - 1)^(delta - 1) * (delta * (p1^-(teta + 1) * 
        teta))))))

                  
 c.copula2.be2 <- (1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * ((p2^-teta - 1)^(delta - 1) * (delta * 
    (p2^-(teta + 1 + 1) * (teta + 1) * teta)) + (p2^-teta - 1)^((delta - 
    1) - 1) * ((delta - 1) * (p2^-(teta + 1) * teta)) * (delta * 
    (p2^-(teta + 1) * teta)))) + ((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((p2^-teta - 
    1)^(delta - 1) * (delta * (p2^-(teta + 1) * teta)))) * ((1/delta) * 
    ((p2^-teta - 1)^(delta - 1) * (delta * (p2^-(teta + 1) * 
        teta)))))) + (1 + ((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^(1/delta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * ((p2^-teta - 1)^(delta - 1) * (delta * 
    (p2^-(teta + 1) * teta)))))) * ((-1/teta) * (((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    ((p2^-teta - 1)^(delta - 1) * (delta * (p2^-(teta + 1) * 
        teta))))))



c.copula2.be1be2 <- (1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((p2^-teta - 1)^(delta - 1) * 
    (delta * (p2^-(teta + 1) * teta)))) * ((1/delta) * ((p1^-teta - 
    1)^(delta - 1) * (delta * (p1^-(teta + 1) * teta)))))) + 
    (1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * (((p1^-teta - 1)^delta + 
        (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
        ((p2^-teta - 1)^(delta - 1) * (delta * (p2^-(teta + 1) * 
            teta)))))) * ((-1/teta) * (((p1^-teta - 1)^delta + 
        (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
        ((p1^-teta - 1)^(delta - 1) * (delta * (p1^-(teta + 1) * 
            teta))))))


c.copula2.be1th <-(-(((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * (log((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))) * 
    (1/teta^2)) - (1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * ((p2^-teta - 1)^(delta - 
    1) * (delta * (p2^-teta * log(p2))) + (p1^-teta - 1)^(delta - 
    1) * (delta * (p1^-teta * log(p1)))))))) * ((-1/teta) * (((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    ((p1^-teta - 1)^(delta - 1) * (delta * (p1^-(teta + 1) * 
        teta)))))) + (1 + ((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^(1/delta))^((-1/teta) - 1) * (1/teta^2 * (((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    ((p1^-teta - 1)^(delta - 1) * (delta * (p1^-(teta + 1) * 
        teta))))) + (-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * ((p1^-teta - 1)^(delta - 
    1) * (delta * (p1^-(teta + 1) - p1^-(teta + 1) * log(p1) * 
    teta)) - (p1^-teta - 1)^((delta - 1) - 1) * ((delta - 1) * 
    (p1^-teta * log(p1))) * (delta * (p1^-(teta + 1) * teta)))) - 
    ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((p2^-teta - 1)^(delta - 
        1) * (delta * (p2^-teta * log(p2))) + (p1^-teta - 1)^(delta - 
        1) * (delta * (p1^-teta * log(p1))))) * ((1/delta) * 
        ((p1^-teta - 1)^(delta - 1) * (delta * (p1^-(teta + 1) * 
            teta)))))))
)*exp(teta.st)



c.copula2.be2th <-(-(((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * (log((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))) * 
    (1/teta^2)) - (1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * ((p2^-teta - 1)^(delta - 
    1) * (delta * (p2^-teta * log(p2))) + (p1^-teta - 1)^(delta - 
    1) * (delta * (p1^-teta * log(p1)))))))) * ((-1/teta) * (((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    ((p2^-teta - 1)^(delta - 1) * (delta * (p2^-(teta + 1) * 
        teta)))))) + (1 + ((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^(1/delta))^((-1/teta) - 1) * (1/teta^2 * (((p1^-teta - 
    1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    ((p2^-teta - 1)^(delta - 1) * (delta * (p2^-(teta + 1) * 
        teta))))) + (-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * ((p2^-teta - 1)^(delta - 
    1) * (delta * (p2^-(teta + 1) - p2^-(teta + 1) * log(p2) * 
    teta)) - (p2^-teta - 1)^((delta - 1) - 1) * ((delta - 1) * 
    (p2^-teta * log(p2))) * (delta * (p2^-(teta + 1) * teta)))) - 
    ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((p2^-teta - 1)^(delta - 
        1) * (delta * (p2^-teta * log(p2))) + (p1^-teta - 1)^(delta - 
        1) * (delta * (p1^-teta * log(p1))))) * ((1/delta) * 
        ((p2^-teta - 1)^(delta - 1) * (delta * (p2^-(teta + 1) * 
            teta))))))))*exp(teta.st)


bit1.th2 <-((1 + ((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 1)^delta)^(1/delta))^(-1/exp(teta.st)) * 
    (log((1 + ((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 
        1)^delta)^(1/delta))) * (exp(teta.st)/exp(teta.st)^2)) - 
    (1 + ((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 
        1)^delta)^(1/delta))^((-1/exp(teta.st)) - 1) * ((-1/exp(teta.st)) * 
        (((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 
            1)^delta)^((1/delta) - 1) * ((1/delta) * ((p2^-exp(teta.st) - 
            1)^(delta - 1) * (delta * (p2^-exp(teta.st) * (log(p2) * 
            exp(teta.st)))) + (p1^-exp(teta.st) - 1)^(delta - 
            1) * (delta * (p1^-exp(teta.st) * (log(p1) * exp(teta.st))))))))) * 
    (log((1 + ((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 
        1)^delta)^(1/delta))) * (exp(teta.st)/exp(teta.st)^2)) + 
    (1 + ((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 
        1)^delta)^(1/delta))^(-1/exp(teta.st)) * (log((1 + ((p1^-exp(teta.st) - 
        1)^delta + (p2^-exp(teta.st) - 1)^delta)^(1/delta))) * 
        (exp(teta.st)/exp(teta.st)^2 - exp(teta.st) * (2 * (exp(teta.st) * 
            exp(teta.st)))/(exp(teta.st)^2)^2) - ((p1^-exp(teta.st) - 
        1)^delta + (p2^-exp(teta.st) - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * ((p2^-exp(teta.st) - 1)^(delta - 1) * 
        (delta * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + 
        (p1^-exp(teta.st) - 1)^(delta - 1) * (delta * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st))))))/(1 + ((p1^-exp(teta.st) - 
        1)^delta + (p2^-exp(teta.st) - 1)^delta)^(1/delta)) * 
        (exp(teta.st)/exp(teta.st)^2)) - (((1 + ((p1^-exp(teta.st) - 
    1)^delta + (p2^-exp(teta.st) - 1)^delta)^(1/delta))^((-1/exp(teta.st)) - 
    1) * (log((1 + ((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 
    1)^delta)^(1/delta))) * (exp(teta.st)/exp(teta.st)^2)) - 
    (1 + ((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 
        1)^delta)^(1/delta))^(((-1/exp(teta.st)) - 1) - 1) * 
        (((-1/exp(teta.st)) - 1) * (((p1^-exp(teta.st) - 1)^delta + 
            (p2^-exp(teta.st) - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
            ((p2^-exp(teta.st) - 1)^(delta - 1) * (delta * (p2^-exp(teta.st) * 
                (log(p2) * exp(teta.st)))) + (p1^-exp(teta.st) - 
                1)^(delta - 1) * (delta * (p1^-exp(teta.st) * 
                (log(p1) * exp(teta.st))))))))) * ((-1/exp(teta.st)) * 
    (((p1^-exp(teta.st) - 1)^delta + (p2^-exp(teta.st) - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * ((p2^-exp(teta.st) - 1)^(delta - 1) * 
        (delta * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + 
        (p1^-exp(teta.st) - 1)^(delta - 1) * (delta * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st)))))))) + (1 + ((p1^-exp(teta.st) - 
    1)^delta + (p2^-exp(teta.st) - 1)^delta)^(1/delta))^((-1/exp(teta.st)) - 
    1) * (exp(teta.st)/exp(teta.st)^2 * (((p1^-exp(teta.st) - 
    1)^delta + (p2^-exp(teta.st) - 1)^delta)^((1/delta) - 1) * 
    ((1/delta) * ((p2^-exp(teta.st) - 1)^(delta - 1) * (delta * 
        (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + (p1^-exp(teta.st) - 
        1)^(delta - 1) * (delta * (p1^-exp(teta.st) * (log(p1) * 
        exp(teta.st))))))) + (-1/exp(teta.st)) * (((p1^-exp(teta.st) - 
    1)^delta + (p2^-exp(teta.st) - 1)^delta)^((1/delta) - 1) * 
    ((1/delta) * ((p2^-exp(teta.st) - 1)^(delta - 1) * (delta * 
        (p2^-exp(teta.st) * (log(p2) * exp(teta.st)) - p2^-exp(teta.st) * 
            (log(p2) * exp(teta.st)) * (log(p2) * exp(teta.st)))) - 
        (p2^-exp(teta.st) - 1)^((delta - 1) - 1) * ((delta - 
            1) * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) * 
            (delta * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + 
        ((p1^-exp(teta.st) - 1)^(delta - 1) * (delta * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st)) - p1^-exp(teta.st) * (log(p1) * 
            exp(teta.st)) * (log(p1) * exp(teta.st)))) - (p1^-exp(teta.st) - 
            1)^((delta - 1) - 1) * ((delta - 1) * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st)))) * (delta * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st))))))) - ((p1^-exp(teta.st) - 
    1)^delta + (p2^-exp(teta.st) - 1)^delta)^(((1/delta) - 1) - 
    1) * (((1/delta) - 1) * ((p2^-exp(teta.st) - 1)^(delta - 
    1) * (delta * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + 
    (p1^-exp(teta.st) - 1)^(delta - 1) * (delta * (p1^-exp(teta.st) * 
        (log(p1) * exp(teta.st)))))) * ((1/delta) * ((p2^-exp(teta.st) - 
    1)^(delta - 1) * (delta * (p2^-exp(teta.st) * (log(p2) * 
    exp(teta.st)))) + (p1^-exp(teta.st) - 1)^(delta - 1) * (delta * 
    (p1^-exp(teta.st) * (log(p1) * exp(teta.st)))))))))

bit1.del2 <-(1 + ((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((p1^-teta - 1)^(exp(delta.st) + 1) * (log((p1^-teta - 
    1)) * exp(delta.st)) + (p2^-teta - 1)^(exp(delta.st) + 1) * 
    (log((p2^-teta - 1)) * exp(delta.st)))) - ((p1^-teta - 1)^(exp(delta.st) + 
    1) + (p2^-teta - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log(((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 
    1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)))) * ((-1/teta) * (((p1^-teta - 1)^(exp(delta.st) + 
    1) + (p2^-teta - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p1^-teta - 1)^(exp(delta.st) + 
    1) * (log((p1^-teta - 1)) * exp(delta.st)) + (p2^-teta - 
    1)^(exp(delta.st) + 1) * (log((p2^-teta - 1)) * exp(delta.st)))) - 
    ((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 1)^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log(((p1^-teta - 1)^(exp(delta.st) + 
        1) + (p2^-teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2)))) + (1 + ((p1^-teta - 1)^(exp(delta.st) + 1) + 
    (p2^-teta - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^((-1/teta) - 
    1) * ((-1/teta) * ((((p1^-teta - 1)^(exp(delta.st) + 1) + 
    (p2^-teta - 1)^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 
    1)) - 1) - 1) * (((1/(exp(delta.st) + 1)) - 1) * ((p1^-teta - 
    1)^(exp(delta.st) + 1) * (log((p1^-teta - 1)) * exp(delta.st)) + 
    (p2^-teta - 1)^(exp(delta.st) + 1) * (log((p2^-teta - 1)) * 
        exp(delta.st)))) - ((p1^-teta - 1)^(exp(delta.st) + 1) + 
    (p2^-teta - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * (log(((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 
    1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * ((1/(exp(delta.st) + 1)) * ((p1^-teta - 1)^(exp(delta.st) + 
    1) * (log((p1^-teta - 1)) * exp(delta.st)) + (p2^-teta - 
    1)^(exp(delta.st) + 1) * (log((p2^-teta - 1)) * exp(delta.st)))) + 
    ((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * ((p1^-teta - 1)^(exp(delta.st) + 1) * (log((p1^-teta - 
        1)) * exp(delta.st)) * (log((p1^-teta - 1)) * exp(delta.st)) + 
        (p1^-teta - 1)^(exp(delta.st) + 1) * (log((p1^-teta - 
            1)) * exp(delta.st)) + ((p2^-teta - 1)^(exp(delta.st) + 
        1) * (log((p2^-teta - 1)) * exp(delta.st)) * (log((p2^-teta - 
        1)) * exp(delta.st)) + (p2^-teta - 1)^(exp(delta.st) + 
        1) * (log((p2^-teta - 1)) * exp(delta.st)))) - exp(delta.st)/(exp(delta.st) + 
        1)^2 * ((p1^-teta - 1)^(exp(delta.st) + 1) * (log((p1^-teta - 
        1)) * exp(delta.st)) + (p2^-teta - 1)^(exp(delta.st) + 
        1) * (log((p2^-teta - 1)) * exp(delta.st)))) - ((((p1^-teta - 
    1)^(exp(delta.st) + 1) + (p2^-teta - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((p1^-teta - 1)^(exp(delta.st) + 1) * (log((p1^-teta - 
    1)) * exp(delta.st)) + (p2^-teta - 1)^(exp(delta.st) + 1) * 
    (log((p2^-teta - 1)) * exp(delta.st)))) - ((p1^-teta - 1)^(exp(delta.st) + 
    1) + (p2^-teta - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log(((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 
    1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * (log(((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 
    1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)) + ((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (((p1^-teta - 
    1)^(exp(delta.st) + 1) * (log((p1^-teta - 1)) * exp(delta.st)) + 
    (p2^-teta - 1)^(exp(delta.st) + 1) * (log((p2^-teta - 1)) * 
        exp(delta.st)))/((p1^-teta - 1)^(exp(delta.st) + 1) + 
    (p2^-teta - 1)^(exp(delta.st) + 1)) * (exp(delta.st)/(exp(delta.st) + 
    1)^2) + log(((p1^-teta - 1)^(exp(delta.st) + 1) + (p2^-teta - 
    1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2 - exp(delta.st) * (2 * (exp(delta.st) * (exp(delta.st) + 
    1)))/((exp(delta.st) + 1)^2)^2)))))


c.copula2.be1del <- (-((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * ((p1^-teta - 1)^delta * 
    log((p1^-teta - 1)) + (p2^-teta - 1)^delta * log((p2^-teta - 
    1)))) - ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta) * 
    (log(((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)) * (1/delta^2)))) * 
    ((-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * ((p1^-teta - 1)^(delta - 1) * (delta * 
        (p1^-(teta + 1) * teta)))))) + (1 + ((p1^-teta - 1)^delta + 
    (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    ((((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((p1^-teta - 1)^delta * 
        log((p1^-teta - 1)) + (p2^-teta - 1)^delta * log((p2^-teta - 
        1)))) - ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
        1) * (log(((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)) * 
        (1/delta^2))) * ((1/delta) * ((p1^-teta - 1)^(delta - 
        1) * (delta * (p1^-(teta + 1) * teta)))) + ((p1^-teta - 
        1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
        ((p1^-teta - 1)^(delta - 1) * log((p1^-teta - 1)) * (delta * 
            (p1^-(teta + 1) * teta)) + (p1^-teta - 1)^(delta - 
            1) * (p1^-(teta + 1) * teta)) - 1/delta^2 * ((p1^-teta - 
        1)^(delta - 1) * (delta * (p1^-(teta + 1) * teta))))))))*exp(delta.st)

c.copula2.be2del <-(-((1 + ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((p1^-teta - 1)^delta + (p2^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * ((p1^-teta - 1)^delta * 
    log((p1^-teta - 1)) + (p2^-teta - 1)^delta * log((p2^-teta - 
    1)))) - ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(1/delta) * 
    (log(((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)) * (1/delta^2)))) * 
    ((-1/teta) * (((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * ((p2^-teta - 1)^(delta - 1) * (delta * 
        (p2^-(teta + 1) * teta)))))) + (1 + ((p1^-teta - 1)^delta + 
    (p2^-teta - 1)^delta)^(1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    ((((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((p1^-teta - 1)^delta * 
        log((p1^-teta - 1)) + (p2^-teta - 1)^delta * log((p2^-teta - 
        1)))) - ((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 
        1) * (log(((p1^-teta - 1)^delta + (p2^-teta - 1)^delta)) * 
        (1/delta^2))) * ((1/delta) * ((p2^-teta - 1)^(delta - 
        1) * (delta * (p2^-(teta + 1) * teta)))) + ((p1^-teta - 
        1)^delta + (p2^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
        ((p2^-teta - 1)^(delta - 1) * log((p2^-teta - 1)) * (delta * 
            (p2^-(teta + 1) * teta)) + (p2^-teta - 1)^(delta - 
            1) * (p2^-(teta + 1) * teta)) - 1/delta^2 * ((p2^-teta - 
        1)^(delta - 1) * (delta * (p2^-(teta + 1) * teta))))))))*exp(delta.st)


bit1.thdel <-(1 + ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^((-1/exp(teta.st)) - 
    1) * ((-1/exp(teta.st)) * (((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) * (log((p1^-exp(teta.st) - 1)) * exp(delta.st)) + 
    (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1) * (log((p2^-exp(teta.st) - 
        1)) * exp(delta.st)))) - ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log(((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + 
    (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)))) * (log((1 + ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))) * (exp(teta.st)/exp(teta.st)^2)) + (1 + ((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^(-1/exp(teta.st)) * ((((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) * (log((p1^-exp(teta.st) - 
    1)) * exp(delta.st)) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) * (log((p2^-exp(teta.st) - 1)) * exp(delta.st)))) - ((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)) * (log(((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)))/(1 + ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + 
    (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (exp(teta.st)/exp(teta.st)^2)) - ((1 + ((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^(((-1/exp(teta.st)) - 1) - 1) * 
    (((-1/exp(teta.st)) - 1) * (((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p1^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) * (log((p1^-exp(teta.st) - 1)) * 
        exp(delta.st)) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) * (log((p2^-exp(teta.st) - 1)) * exp(delta.st)))) - 
        ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
            1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
            (log(((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + 
                (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * 
                (exp(delta.st)/(exp(delta.st) + 1)^2)))) * ((-1/exp(teta.st)) * 
    (((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
        1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 1) * 
        ((1/(exp(delta.st) + 1)) * ((p2^-exp(teta.st) - 1)^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * (p2^-exp(teta.st) * 
            (log(p2) * exp(teta.st)))) + (p1^-exp(teta.st) - 
            1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
            1) * (p1^-exp(teta.st) * (log(p1) * exp(teta.st)))))))) + 
    (1 + ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
        1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^((-1/exp(teta.st)) - 
        1) * ((-1/exp(teta.st)) * ((((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 
        1)) - 1) - 1) * (((1/(exp(delta.st) + 1)) - 1) * ((p1^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) * (log((p1^-exp(teta.st) - 1)) * 
        exp(delta.st)) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) * (log((p2^-exp(teta.st) - 1)) * exp(delta.st)))) - 
        ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
            1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
            1) * (log(((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
            1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * 
            (exp(delta.st)/(exp(delta.st) + 1)^2))) * ((1/(exp(delta.st) + 
        1)) * ((p2^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 
        1) * ((exp(delta.st) + 1) * (p2^-exp(teta.st) * (log(p2) * 
        exp(teta.st)))) + (p1^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * ((exp(delta.st) + 1) * (p1^-exp(teta.st) * 
        (log(p1) * exp(teta.st)))))) + ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p2^-exp(teta.st) - 
        1)^((exp(delta.st) + 1) - 1) * (log((p2^-exp(teta.st) - 
        1)) * exp(delta.st)) * ((exp(delta.st) + 1) * (p2^-exp(teta.st) * 
        (log(p2) * exp(teta.st)))) + (p2^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * (exp(delta.st) * (p2^-exp(teta.st) * (log(p2) * 
        exp(teta.st)))) + ((p1^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * (log((p1^-exp(teta.st) - 1)) * exp(delta.st)) * 
        ((exp(delta.st) + 1) * (p1^-exp(teta.st) * (log(p1) * 
            exp(teta.st)))) + (p1^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * (exp(delta.st) * (p1^-exp(teta.st) * (log(p1) * 
        exp(teta.st)))))) - exp(delta.st)/(exp(delta.st) + 1)^2 * 
        ((p2^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
            1) * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + 
            (p1^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 1) * 
                ((exp(delta.st) + 1) * (p1^-exp(teta.st) * (log(p1) * 
                  exp(teta.st)))))))))


}


if(BivD=="BB1.90"){


  c.copula.be1 <- -((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - p1)^teta - 
    1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 1) * teta)))))))


 
  c.copula.be2 <- 1 + (1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 
    1) * (delta * (p2^(teta - 1) * teta))))))

  c.copula.theta <- ((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^(1/teta) * 
    (log((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))) * 
        (1/teta^2)) + (1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^(-1/delta))^((1/teta) - 1) * ((1/teta) * ((((1 - 
    p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 1) * (delta * 
    (p2^teta * log(p2))) + ((1 - p1)^teta - 1)^-(delta + 1) * 
    (delta * ((1 - p1)^teta * log((1 - p1)))))))))*(-exp(teta.st))


  c.copula.delta <- (-((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^(-1/delta) * (log((((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)) * (1/delta^2)) - (((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * 
    ((p2^teta - 1)^-delta * log((p2^teta - 1)) + ((1 - p1)^teta - 
        1)^-delta * log(((1 - p1)^teta - 1))))))))*(-exp(delta.st))


 c.copula2.be1 <-  -((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - 
    p1)^teta - 1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 1) * 
    teta)))))) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - p1)^teta - 
    1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 1) * teta)))))) + 
    (1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
        1) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
        1)^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * 
        (((1 - p1)^teta - 1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
            1) * teta)))) * ((-1/delta) * (((1 - p1)^teta - 1)^-(delta + 
        1) * (delta * ((1 - p1)^(teta - 1) * teta)))) + (((1 - 
        p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * (((1 - p1)^teta - 1)^-(delta + 1 + 
        1) * ((delta + 1) * ((1 - p1)^(teta - 1) * teta)) * (delta * 
        ((1 - p1)^(teta - 1) * teta)) - ((1 - p1)^teta - 1)^-(delta + 
        1) * (delta * ((1 - p1)^((teta - 1) - 1) * (teta - 1) * 
        teta)))))))

                  
 c.copula2.be2 <- (1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 
    1) * (delta * (p2^((teta - 1) - 1) * (teta - 1) * teta)) - 
    (p2^teta - 1)^-(delta + 1 + 1) * ((delta + 1) * (p2^(teta - 
        1) * teta)) * (delta * (p2^(teta - 1) * teta)))) - (((1 - 
    p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(((-1/delta) - 
    1) - 1) * (((-1/delta) - 1) * ((p2^teta - 1)^-(delta + 1) * 
    (delta * (p2^(teta - 1) * teta)))) * ((-1/delta) * ((p2^teta - 
    1)^-(delta + 1) * (delta * (p2^(teta - 1) * teta)))))) - 
    (1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * ((((1 - p1)^teta - 1)^-delta + 
        (p2^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((p2^teta - 1)^-(delta + 1) * (delta * (p2^(teta - 1) * 
            teta)))))) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + 
        (p2^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((p2^teta - 1)^-(delta + 1) * (delta * (p2^(teta - 1) * 
            teta))))))


                  
                 

c.copula2.be1be2 <- (1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((p2^teta - 
    1)^-(delta + 1) * (delta * (p2^(teta - 1) * teta)))) * ((-1/delta) * 
    (((1 - p1)^teta - 1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
        1) * teta)))))) + (1 + (((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 1) - 1) * 
    (((1/teta) - 1) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
        1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p2^teta - 
        1)^-(delta + 1) * (delta * (p2^(teta - 1) * teta)))))) * 
    ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * (((1 - p1)^teta - 1)^-(delta + 1) * 
        (delta * ((1 - p1)^(teta - 1) * teta))))))



c.copula2.be1th <-(-((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - p1)^teta - 
    1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 1) * log((1 - 
    p1)) * teta + (1 - p1)^(teta - 1))) - ((1 - p1)^teta - 1)^-(delta + 
    1 + 1) * ((delta + 1) * ((1 - p1)^teta * log((1 - p1)))) * 
    (delta * ((1 - p1)^(teta - 1) * teta)))) - (((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)^(((-1/delta) - 1) - 1) * 
    (((-1/delta) - 1) * ((p2^teta - 1)^-(delta + 1) * (delta * 
        (p2^teta * log(p2))) + ((1 - p1)^teta - 1)^-(delta + 
        1) * (delta * ((1 - p1)^teta * log((1 - p1)))))) * ((-1/delta) * 
    (((1 - p1)^teta - 1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
        1) * teta))))) - 1/teta^2 * ((((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - 
    p1)^teta - 1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 1) * 
    teta)))))) - ((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^(-1/delta))^((1/teta) - 1) * (log((1 + (((1 - 
    p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))) * 
    (1/teta^2)) + (1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 
    1) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 1) * (delta * 
    (p2^teta * log(p2))) + ((1 - p1)^teta - 1)^-(delta + 1) * 
    (delta * ((1 - p1)^teta * log((1 - p1))))))))) * ((1/teta) * 
    ((((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * (((1 - p1)^teta - 1)^-(delta + 1) * 
        (delta * ((1 - p1)^(teta - 1) * teta))))))))*(-exp(teta.st))


        


c.copula2.be2th <-((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 
    1) * (delta * (p2^(teta - 1) * log(p2) * teta + p2^(teta - 
    1))) - (p2^teta - 1)^-(delta + 1 + 1) * ((delta + 1) * (p2^teta * 
    log(p2))) * (delta * (p2^(teta - 1) * teta)))) - (((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)^(((-1/delta) - 1) - 1) * 
    (((-1/delta) - 1) * ((p2^teta - 1)^-(delta + 1) * (delta * 
        (p2^teta * log(p2))) + ((1 - p1)^teta - 1)^-(delta + 
        1) * (delta * ((1 - p1)^teta * log((1 - p1)))))) * ((-1/delta) * 
    ((p2^teta - 1)^-(delta + 1) * (delta * (p2^(teta - 1) * teta))))) - 
    1/teta^2 * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 1) * (delta * 
        (p2^(teta - 1) * teta)))))) - ((1 + (((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * (log((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))) * 
    (1/teta^2)) + (1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 
    1) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 1) * (delta * 
    (p2^teta * log(p2))) + ((1 - p1)^teta - 1)^-(delta + 1) * 
    (delta * ((1 - p1)^teta * log((1 - p1))))))))) * ((1/teta) * 
    ((((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 1) * (delta * 
        (p2^(teta - 1) * teta)))))))*(-exp(teta.st))






bit1.th2 <--(((1 + (((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
    1)^-delta)^(-1/delta))^(((-1/exp(teta.st)) - 1) - 1) * (((-1/exp(teta.st)) - 
    1) * ((((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - p1)^-exp(teta.st) - 
    1)^-(delta + 1) * (delta * ((1 - p1)^-exp(teta.st) * (log((1 - 
    p1)) * exp(teta.st)))) + (p2^-exp(teta.st) - 1)^-(delta + 
    1) * (delta * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))))))) + 
    (1 + (((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
        1)^-delta)^(-1/delta))^((-1/exp(teta.st)) - 1) * (log((1 + 
        (((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
            1)^-delta)^(-1/delta))) * (exp(teta.st)/exp(teta.st)^2))) * 
    ((-1/exp(teta.st)) * ((((1 - p1)^-exp(teta.st) - 1)^-delta + 
        (p2^-exp(teta.st) - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        (((1 - p1)^-exp(teta.st) - 1)^-(delta + 1) * (delta * 
            ((1 - p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st)))) + 
            (p2^-exp(teta.st) - 1)^-(delta + 1) * (delta * (p2^-exp(teta.st) * 
                (log(p2) * exp(teta.st)))))))) + (1 + (((1 - 
    p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 1)^-delta)^(-1/delta))^((-1/exp(teta.st)) - 
    1) * (exp(teta.st)/exp(teta.st)^2 * ((((1 - p1)^-exp(teta.st) - 
    1)^-delta + (p2^-exp(teta.st) - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * (((1 - p1)^-exp(teta.st) - 1)^-(delta + 
    1) * (delta * ((1 - p1)^-exp(teta.st) * (log((1 - p1)) * 
    exp(teta.st)))) + (p2^-exp(teta.st) - 1)^-(delta + 1) * (delta * 
    (p2^-exp(teta.st) * (log(p2) * exp(teta.st))))))) + (-1/exp(teta.st)) * 
    ((((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
        1)^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * 
        (((1 - p1)^-exp(teta.st) - 1)^-(delta + 1) * (delta * 
            ((1 - p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st)))) + 
            (p2^-exp(teta.st) - 1)^-(delta + 1) * (delta * (p2^-exp(teta.st) * 
                (log(p2) * exp(teta.st)))))) * ((-1/delta) * 
        (((1 - p1)^-exp(teta.st) - 1)^-(delta + 1) * (delta * 
            ((1 - p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st)))) + 
            (p2^-exp(teta.st) - 1)^-(delta + 1) * (delta * (p2^-exp(teta.st) * 
                (log(p2) * exp(teta.st)))))) + (((1 - p1)^-exp(teta.st) - 
        1)^-delta + (p2^-exp(teta.st) - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * (((1 - p1)^-exp(teta.st) - 1)^-(delta + 
        1 + 1) * ((delta + 1) * ((1 - p1)^-exp(teta.st) * (log((1 - 
        p1)) * exp(teta.st)))) * (delta * ((1 - p1)^-exp(teta.st) * 
        (log((1 - p1)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 
        1)^-(delta + 1) * (delta * ((1 - p1)^-exp(teta.st) * 
        (log((1 - p1)) * exp(teta.st)) - (1 - p1)^-exp(teta.st) * 
        (log((1 - p1)) * exp(teta.st)) * (log((1 - p1)) * exp(teta.st)))) + 
        ((p2^-exp(teta.st) - 1)^-(delta + 1 + 1) * ((delta + 
            1) * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) * 
            (delta * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + 
            (p2^-exp(teta.st) - 1)^-(delta + 1) * (delta * (p2^-exp(teta.st) * 
                (log(p2) * exp(teta.st)) - p2^-exp(teta.st) * 
                (log(p2) * exp(teta.st)) * (log(p2) * exp(teta.st))))))))) + 
    (((1 + (((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
        1)^-delta)^(-1/delta))^((-1/exp(teta.st)) - 1) * ((-1/exp(teta.st)) * 
        ((((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
            1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - 
            p1)^-exp(teta.st) - 1)^-(delta + 1) * (delta * ((1 - 
            p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st)))) + 
            (p2^-exp(teta.st) - 1)^-(delta + 1) * (delta * (p2^-exp(teta.st) * 
                (log(p2) * exp(teta.st)))))))) + (1 + (((1 - 
        p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 1)^-delta)^(-1/delta))^(-1/exp(teta.st)) * 
        (log((1 + (((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
            1)^-delta)^(-1/delta))) * (exp(teta.st)/exp(teta.st)^2))) * 
        (log((1 + (((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
            1)^-delta)^(-1/delta))) * (exp(teta.st)/exp(teta.st)^2)) + 
        (1 + (((1 - p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
            1)^-delta)^(-1/delta))^(-1/exp(teta.st)) * ((((1 - 
            p1)^-exp(teta.st) - 1)^-delta + (p2^-exp(teta.st) - 
            1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - 
            p1)^-exp(teta.st) - 1)^-(delta + 1) * (delta * ((1 - 
            p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st)))) + 
            (p2^-exp(teta.st) - 1)^-(delta + 1) * (delta * (p2^-exp(teta.st) * 
                (log(p2) * exp(teta.st))))))/(1 + (((1 - p1)^-exp(teta.st) - 
            1)^-delta + (p2^-exp(teta.st) - 1)^-delta)^(-1/delta)) * 
            (exp(teta.st)/exp(teta.st)^2) + log((1 + (((1 - p1)^-exp(teta.st) - 
            1)^-delta + (p2^-exp(teta.st) - 1)^-delta)^(-1/delta))) * 
            (exp(teta.st)/exp(teta.st)^2 - exp(teta.st) * (2 * 
                (exp(teta.st) * exp(teta.st)))/(exp(teta.st)^2)^2))))





bit1.del2 <--((1 + (((1 - p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((((1 - p1)^teta - 1)^(exp(delta.st) + 
    1) + (p2^teta - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * (((1 - p1)^teta - 1)^(exp(delta.st) + 
    1) * (log(((1 - p1)^teta - 1)) * exp(delta.st)) + (p2^teta - 
    1)^(exp(delta.st) + 1) * (log((p2^teta - 1)) * exp(delta.st)))) - 
    (((1 - p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log((((1 - p1)^teta - 
        1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)))) * ((1/teta) * 
    ((((1 - p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * (((1 - p1)^teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
        p1)^teta - 1)) * exp(delta.st)) + (p2^teta - 1)^(exp(delta.st) + 
        1) * (log((p2^teta - 1)) * exp(delta.st)))) - (((1 - 
        p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log((((1 - p1)^teta - 
        1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)))) + (1 + 
    (((1 - p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)))^((1/teta) - 1) * ((1/teta) * 
    (((((1 - p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^(((1/(exp(delta.st) + 1)) - 1) - 1) * (((1/(exp(delta.st) + 
        1)) - 1) * (((1 - p1)^teta - 1)^(exp(delta.st) + 1) * 
        (log(((1 - p1)^teta - 1)) * exp(delta.st)) + (p2^teta - 
        1)^(exp(delta.st) + 1) * (log((p2^teta - 1)) * exp(delta.st)))) - 
        (((1 - p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 
            1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
            1) * (log((((1 - p1)^teta - 1)^(exp(delta.st) + 1) + 
            (p2^teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
            1)^2))) * ((1/(exp(delta.st) + 1)) * (((1 - p1)^teta - 
        1)^(exp(delta.st) + 1) * (log(((1 - p1)^teta - 1)) * 
        exp(delta.st)) + (p2^teta - 1)^(exp(delta.st) + 1) * 
        (log((p2^teta - 1)) * exp(delta.st)))) + (((1 - p1)^teta - 
        1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * (((1 - p1)^teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
        p1)^teta - 1)) * exp(delta.st)) * (log(((1 - p1)^teta - 
        1)) * exp(delta.st)) + ((1 - p1)^teta - 1)^(exp(delta.st) + 
        1) * (log(((1 - p1)^teta - 1)) * exp(delta.st)) + ((p2^teta - 
        1)^(exp(delta.st) + 1) * (log((p2^teta - 1)) * exp(delta.st)) * 
        (log((p2^teta - 1)) * exp(delta.st)) + (p2^teta - 1)^(exp(delta.st) + 
        1) * (log((p2^teta - 1)) * exp(delta.st)))) - exp(delta.st)/(exp(delta.st) + 
        1)^2 * (((1 - p1)^teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
        p1)^teta - 1)) * exp(delta.st)) + (p2^teta - 1)^(exp(delta.st) + 
        1) * (log((p2^teta - 1)) * exp(delta.st)))) - (((((1 - 
        p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * (((1 - p1)^teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
        p1)^teta - 1)) * exp(delta.st)) + (p2^teta - 1)^(exp(delta.st) + 
        1) * (log((p2^teta - 1)) * exp(delta.st)))) - (((1 - 
        p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log((((1 - p1)^teta - 
        1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * (log((((1 - 
        p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)) + (((1 - 
        p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * ((((1 - p1)^teta - 1)^(exp(delta.st) + 
        1) * (log(((1 - p1)^teta - 1)) * exp(delta.st)) + (p2^teta - 
        1)^(exp(delta.st) + 1) * (log((p2^teta - 1)) * exp(delta.st)))/(((1 - 
        p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1)) * (exp(delta.st)/(exp(delta.st) + 1)^2) + log((((1 - 
        p1)^teta - 1)^(exp(delta.st) + 1) + (p2^teta - 1)^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2 - exp(delta.st) * 
        (2 * (exp(delta.st) * (exp(delta.st) + 1)))/((exp(delta.st) + 
        1)^2)^2))))))




c.copula2.be1del <- (-((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)^(-1/delta) * (log((((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)) * (1/delta^2)) - (((1 - 
    p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((p2^teta - 1)^-delta * log((p2^teta - 
    1)) + ((1 - p1)^teta - 1)^-delta * log(((1 - p1)^teta - 1)))))) * 
    ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * (((1 - p1)^teta - 1)^-(delta + 1) * 
        (delta * ((1 - p1)^(teta - 1) * teta)))))) + (1 + (((1 - 
    p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^((-1/delta) - 1) * (log((((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)) * (1/delta^2)) - (((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)^(((-1/delta) - 1) - 1) * 
    (((-1/delta) - 1) * ((p2^teta - 1)^-delta * log((p2^teta - 
        1)) + ((1 - p1)^teta - 1)^-delta * log(((1 - p1)^teta - 
        1))))) * ((-1/delta) * (((1 - p1)^teta - 1)^-(delta + 
    1) * (delta * ((1 - p1)^(teta - 1) * teta)))) + (((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 1) * (1/delta^2 * 
    (((1 - p1)^teta - 1)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
        1) * teta))) + (-1/delta) * (((1 - p1)^teta - 1)^-(delta + 
    1) * ((1 - p1)^(teta - 1) * teta) - ((1 - p1)^teta - 1)^-(delta + 
    1) * log(((1 - p1)^teta - 1)) * (delta * ((1 - p1)^(teta - 
    1) * teta))))))))*(-exp(delta.st))


        



        
c.copula2.be2del <- ((1 + (((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)^(-1/delta) * (log((((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)) * (1/delta^2)) - (((1 - 
    p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((p2^teta - 1)^-delta * log((p2^teta - 
    1)) + ((1 - p1)^teta - 1)^-delta * log(((1 - p1)^teta - 1)))))) * 
    ((1/teta) * ((((1 - p1)^teta - 1)^-delta + (p2^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((p2^teta - 1)^-(delta + 1) * (delta * 
        (p2^(teta - 1) * teta)))))) + (1 + (((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((((1 - p1)^teta - 1)^-delta + (p2^teta - 
    1)^-delta)^((-1/delta) - 1) * (log((((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)) * (1/delta^2)) - (((1 - p1)^teta - 
    1)^-delta + (p2^teta - 1)^-delta)^(((-1/delta) - 1) - 1) * 
    (((-1/delta) - 1) * ((p2^teta - 1)^-delta * log((p2^teta - 
        1)) + ((1 - p1)^teta - 1)^-delta * log(((1 - p1)^teta - 
        1))))) * ((-1/delta) * ((p2^teta - 1)^-(delta + 1) * 
    (delta * (p2^(teta - 1) * teta)))) + (((1 - p1)^teta - 1)^-delta + 
    (p2^teta - 1)^-delta)^((-1/delta) - 1) * (1/delta^2 * ((p2^teta - 
    1)^-(delta + 1) * (delta * (p2^(teta - 1) * teta))) + (-1/delta) * 
    ((p2^teta - 1)^-(delta + 1) * (p2^(teta - 1) * teta) - (p2^teta - 
        1)^-(delta + 1) * log((p2^teta - 1)) * (delta * (p2^(teta - 
        1) * teta)))))))*(-exp(delta.st))




bit1.thdel <--((1 + (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^((1/-exp(teta.st)) - 
    1) * ((1/-exp(teta.st)) * ((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * (((1 - p1)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) * (log(((1 - p1)^-exp(teta.st) - 1)) * 
    exp(delta.st)) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) * (log((p2^-exp(teta.st) - 1)) * exp(delta.st)))) - (((1 - 
    p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log((((1 - 
    p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)))) * (log((1 + (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))) * (exp(teta.st)/(-exp(teta.st))^2)) + (1 + (((1 - p1)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^(1/-exp(teta.st)) * (((((1 - 
    p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) * 
    (log(((1 - p1)^-exp(teta.st) - 1)) * exp(delta.st)) + (p2^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) * (log((p2^-exp(teta.st) - 1)) * exp(delta.st)))) - 
    (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
        1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log((((1 - 
        p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
        1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2)))/(1 + (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (exp(teta.st)/(-exp(teta.st))^2)) - ((1 + (((1 - p1)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^(((1/-exp(teta.st)) - 1) - 1) * 
    (((1/-exp(teta.st)) - 1) * ((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * (((1 - p1)^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) * (log(((1 - p1)^-exp(teta.st) - 
        1)) * exp(delta.st)) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) * (log((p2^-exp(teta.st) - 1)) * exp(delta.st)))) - 
        (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
            1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
            (log((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
                1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)))) * 
    ((1/-exp(teta.st)) * ((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p2^-exp(teta.st) - 
        1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * 
        (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + ((1 - 
        p1)^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
        1) * ((1 - p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st)))))))) + 
    (1 + (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + 
        (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)))^((1/-exp(teta.st)) - 1) * ((1/-exp(teta.st)) * (((((1 - 
        p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
        1)^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 1)) - 1) - 
        1) * (((1/(exp(delta.st) + 1)) - 1) * (((1 - p1)^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) * (log(((1 - p1)^-exp(teta.st) - 
        1)) * exp(delta.st)) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) * (log((p2^-exp(teta.st) - 1)) * exp(delta.st)))) - 
        (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 
            1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
            1) * (log((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
            1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * 
            (exp(delta.st)/(exp(delta.st) + 1)^2))) * ((1/(exp(delta.st) + 
        1)) * ((p2^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 
        1) * ((exp(delta.st) + 1) * (p2^-exp(teta.st) * (log(p2) * 
        exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * ((exp(delta.st) + 1) * ((1 - p1)^-exp(teta.st) * 
        (log((1 - p1)) * exp(teta.st)))))) + (((1 - p1)^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) + (p2^-exp(teta.st) - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * ((p2^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 
        1) * (log((p2^-exp(teta.st) - 1)) * exp(delta.st)) * 
        ((exp(delta.st) + 1) * (p2^-exp(teta.st) * (log(p2) * 
            exp(teta.st)))) + (p2^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * (exp(delta.st) * (p2^-exp(teta.st) * (log(p2) * 
        exp(teta.st)))) + (((1 - p1)^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * (log(((1 - p1)^-exp(teta.st) - 1)) * exp(delta.st)) * 
        ((exp(delta.st) + 1) * ((1 - p1)^-exp(teta.st) * (log((1 - 
            p1)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 
        1)^((exp(delta.st) + 1) - 1) * (exp(delta.st) * ((1 - 
        p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st)))))) - 
        exp(delta.st)/(exp(delta.st) + 1)^2 * ((p2^-exp(teta.st) - 
            1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
            1) * (p2^-exp(teta.st) * (log(p2) * exp(teta.st)))) + 
            ((1 - p1)^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 
                1) * ((exp(delta.st) + 1) * ((1 - p1)^-exp(teta.st) * 
                (log((1 - p1)) * exp(teta.st))))))))))



}


if(BivD=="BB1.180"){
   
  c.copula.be1 <-1 + (1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p1)^-teta - 
    1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 1) * teta))))))

 
  c.copula.be2 <-1 + (1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p2)^-teta - 
    1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 1) * teta))))))

  c.copula.theta <- ((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^(-1/teta) * 
    (log((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
        1)^delta)^(1/delta))) * (1/teta^2)) - (1 + (((1 - p1)^-teta - 
    1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p2)^-teta - 
    1)^(delta - 1) * (delta * ((1 - p2)^-teta * log((1 - p2)))) + 
    ((1 - p1)^-teta - 1)^(delta - 1) * (delta * ((1 - p1)^-teta * 
        log((1 - p1)))))))))*exp(teta.st)


  c.copula.delta <- ((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p1)^-teta - 
    1)^delta * log(((1 - p1)^-teta - 1)) + ((1 - p2)^-teta - 
    1)^delta * log(((1 - p2)^-teta - 1)))) - (((1 - p1)^-teta - 
    1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta) * (log((((1 - 
    p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)) * (1/delta^2)))))*exp(delta.st)



 c.copula2.be1 <- (1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p1)^-teta - 1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta)))))) * ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p1)^-teta - 1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta)))))) + (1 + (((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 1) * 
    ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
        1)^delta)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * 
        (((1 - p1)^-teta - 1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
            1) * teta)))) * ((1/delta) * (((1 - p1)^-teta - 1)^(delta - 
        1) * (delta * ((1 - p1)^-(teta + 1) * teta)))) + (((1 - 
        p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * (((1 - p1)^-teta - 1)^((delta - 1) - 
        1) * ((delta - 1) * ((1 - p1)^-(teta + 1) * teta)) * 
        (delta * ((1 - p1)^-(teta + 1) * teta)) + ((1 - p1)^-teta - 
        1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 1 + 1) * 
        (teta + 1) * teta))))))


                  
 c.copula2.be2 <-(1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p2)^-teta - 1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta)))))) * ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p2)^-teta - 1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta)))))) + (1 + (((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 1) * 
    ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
        1)^delta)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * 
        (((1 - p2)^-teta - 1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
            1) * teta)))) * ((1/delta) * (((1 - p2)^-teta - 1)^(delta - 
        1) * (delta * ((1 - p2)^-(teta + 1) * teta)))) + (((1 - 
        p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * (((1 - p2)^-teta - 1)^((delta - 1) - 
        1) * ((delta - 1) * ((1 - p2)^-(teta + 1) * teta)) * 
        (delta * ((1 - p2)^-(teta + 1) * teta)) + ((1 - p2)^-teta - 
        1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 1 + 1) * 
        (teta + 1) * teta))))))

                  
                 



c.copula2.be1be2 <- (1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p2)^-teta - 1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta)))))) * ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p1)^-teta - 1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta)))))) + (1 + (((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 1) * 
    ((-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
        1)^delta)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * 
        (((1 - p2)^-teta - 1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
            1) * teta)))) * ((1/delta) * (((1 - p1)^-teta - 1)^(delta - 
        1) * (delta * ((1 - p1)^-(teta + 1) * teta))))))




c.copula2.be1th <-(((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * (log((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^(1/delta))) * (1/teta^2)) - (1 + (((1 - p1)^-teta - 
    1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p2)^-teta - 1)^(delta - 1) * (delta * ((1 - p2)^-teta * 
        log((1 - p2)))) + ((1 - p1)^-teta - 1)^(delta - 1) * 
        (delta * ((1 - p1)^-teta * log((1 - p1))))))))) * ((-1/teta) * 
    ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * (((1 - p1)^-teta - 1)^(delta - 1) * 
        (delta * ((1 - p1)^-(teta + 1) * teta)))))) + (1 + (((1 - 
    p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * (1/teta^2 * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p1)^-teta - 
    1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 1) * teta))))) + 
    (-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
        1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p1)^-teta - 
        1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 1) - (1 - 
        p1)^-(teta + 1) * log((1 - p1)) * teta)) - ((1 - p1)^-teta - 
        1)^((delta - 1) - 1) * ((delta - 1) * ((1 - p1)^-teta * 
        log((1 - p1)))) * (delta * ((1 - p1)^-(teta + 1) * teta)))) - 
        (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(((1/delta) - 
            1) - 1) * (((1/delta) - 1) * (((1 - p2)^-teta - 1)^(delta - 
            1) * (delta * ((1 - p2)^-teta * log((1 - p2)))) + 
            ((1 - p1)^-teta - 1)^(delta - 1) * (delta * ((1 - 
                p1)^-teta * log((1 - p1)))))) * ((1/delta) * 
            (((1 - p1)^-teta - 1)^(delta - 1) * (delta * ((1 - 
                p1)^-(teta + 1) * teta)))))))*exp(teta.st)


        

c.copula2.be2th <-(((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * (log((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^(1/delta))) * (1/teta^2)) - (1 + (((1 - p1)^-teta - 
    1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p2)^-teta - 1)^(delta - 1) * (delta * ((1 - p2)^-teta * 
        log((1 - p2)))) + ((1 - p1)^-teta - 1)^(delta - 1) * 
        (delta * ((1 - p1)^-teta * log((1 - p1))))))))) * ((-1/teta) * 
    ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * (((1 - p2)^-teta - 1)^(delta - 1) * 
        (delta * ((1 - p2)^-(teta + 1) * teta)))))) + (1 + (((1 - 
    p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * (1/teta^2 * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p2)^-teta - 
    1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 1) * teta))))) + 
    (-1/teta) * ((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
        1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p2)^-teta - 
        1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 1) - (1 - 
        p2)^-(teta + 1) * log((1 - p2)) * teta)) - ((1 - p2)^-teta - 
        1)^((delta - 1) - 1) * ((delta - 1) * ((1 - p2)^-teta * 
        log((1 - p2)))) * (delta * ((1 - p2)^-(teta + 1) * teta)))) - 
        (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(((1/delta) - 
            1) - 1) * (((1/delta) - 1) * (((1 - p2)^-teta - 1)^(delta - 
            1) * (delta * ((1 - p2)^-teta * log((1 - p2)))) + 
            ((1 - p1)^-teta - 1)^(delta - 1) * (delta * ((1 - 
                p1)^-teta * log((1 - p1)))))) * ((1/delta) * 
            (((1 - p2)^-teta - 1)^(delta - 1) * (delta * ((1 - 
                p2)^-(teta + 1) * teta)))))))*exp(teta.st)





bit1.th2 <-((1 + (((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
    1)^delta)^(1/delta))^(-1/exp(teta.st)) * (log((1 + (((1 - 
    p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
    1)^delta)^(1/delta))) * (exp(teta.st)/exp(teta.st)^2)) - 
    (1 + (((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
        1)^delta)^(1/delta))^((-1/exp(teta.st)) - 1) * ((-1/exp(teta.st)) * 
        ((((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
            1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p2)^-exp(teta.st) - 
            1)^(delta - 1) * (delta * ((1 - p2)^-exp(teta.st) * 
            (log((1 - p2)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 
            1)^(delta - 1) * (delta * ((1 - p1)^-exp(teta.st) * 
            (log((1 - p1)) * exp(teta.st))))))))) * (log((1 + 
    (((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
        1)^delta)^(1/delta))) * (exp(teta.st)/exp(teta.st)^2)) + 
    (1 + (((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
        1)^delta)^(1/delta))^(-1/exp(teta.st)) * (log((1 + (((1 - 
        p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
        1)^delta)^(1/delta))) * (exp(teta.st)/exp(teta.st)^2 - 
        exp(teta.st) * (2 * (exp(teta.st) * exp(teta.st)))/(exp(teta.st)^2)^2) - 
        (((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
            1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p2)^-exp(teta.st) - 
            1)^(delta - 1) * (delta * ((1 - p2)^-exp(teta.st) * 
            (log((1 - p2)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 
            1)^(delta - 1) * (delta * ((1 - p1)^-exp(teta.st) * 
            (log((1 - p1)) * exp(teta.st))))))/(1 + (((1 - p1)^-exp(teta.st) - 
            1)^delta + ((1 - p2)^-exp(teta.st) - 1)^delta)^(1/delta)) * 
            (exp(teta.st)/exp(teta.st)^2)) - (((1 + (((1 - p1)^-exp(teta.st) - 
    1)^delta + ((1 - p2)^-exp(teta.st) - 1)^delta)^(1/delta))^((-1/exp(teta.st)) - 
    1) * (log((1 + (((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - 
    p2)^-exp(teta.st) - 1)^delta)^(1/delta))) * (exp(teta.st)/exp(teta.st)^2)) - 
    (1 + (((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
        1)^delta)^(1/delta))^(((-1/exp(teta.st)) - 1) - 1) * 
        (((-1/exp(teta.st)) - 1) * ((((1 - p1)^-exp(teta.st) - 
            1)^delta + ((1 - p2)^-exp(teta.st) - 1)^delta)^((1/delta) - 
            1) * ((1/delta) * (((1 - p2)^-exp(teta.st) - 1)^(delta - 
            1) * (delta * ((1 - p2)^-exp(teta.st) * (log((1 - 
            p2)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 
            1)^(delta - 1) * (delta * ((1 - p1)^-exp(teta.st) * 
            (log((1 - p1)) * exp(teta.st))))))))) * ((-1/exp(teta.st)) * 
    ((((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
        1)^delta)^((1/delta) - 1) * ((1/delta) * (((1 - p2)^-exp(teta.st) - 
        1)^(delta - 1) * (delta * ((1 - p2)^-exp(teta.st) * (log((1 - 
        p2)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 1)^(delta - 
        1) * (delta * ((1 - p1)^-exp(teta.st) * (log((1 - p1)) * 
        exp(teta.st)))))))) + (1 + (((1 - p1)^-exp(teta.st) - 
    1)^delta + ((1 - p2)^-exp(teta.st) - 1)^delta)^(1/delta))^((-1/exp(teta.st)) - 
    1) * (exp(teta.st)/exp(teta.st)^2 * ((((1 - p1)^-exp(teta.st) - 
    1)^delta + ((1 - p2)^-exp(teta.st) - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * (((1 - p2)^-exp(teta.st) - 1)^(delta - 
    1) * (delta * ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * 
    exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 1)^(delta - 1) * 
    (delta * ((1 - p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st))))))) + 
    (-1/exp(teta.st)) * ((((1 - p1)^-exp(teta.st) - 1)^delta + 
        ((1 - p2)^-exp(teta.st) - 1)^delta)^((1/delta) - 1) * 
        ((1/delta) * (((1 - p2)^-exp(teta.st) - 1)^(delta - 1) * 
            (delta * ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * 
                exp(teta.st)) - (1 - p2)^-exp(teta.st) * (log((1 - 
                p2)) * exp(teta.st)) * (log((1 - p2)) * exp(teta.st)))) - 
            ((1 - p2)^-exp(teta.st) - 1)^((delta - 1) - 1) * 
                ((delta - 1) * ((1 - p2)^-exp(teta.st) * (log((1 - 
                  p2)) * exp(teta.st)))) * (delta * ((1 - p2)^-exp(teta.st) * 
                (log((1 - p2)) * exp(teta.st)))) + (((1 - p1)^-exp(teta.st) - 
            1)^(delta - 1) * (delta * ((1 - p1)^-exp(teta.st) * 
            (log((1 - p1)) * exp(teta.st)) - (1 - p1)^-exp(teta.st) * 
            (log((1 - p1)) * exp(teta.st)) * (log((1 - p1)) * 
            exp(teta.st)))) - ((1 - p1)^-exp(teta.st) - 1)^((delta - 
            1) - 1) * ((delta - 1) * ((1 - p1)^-exp(teta.st) * 
            (log((1 - p1)) * exp(teta.st)))) * (delta * ((1 - 
            p1)^-exp(teta.st) * (log((1 - p1)) * exp(teta.st))))))) - 
        (((1 - p1)^-exp(teta.st) - 1)^delta + ((1 - p2)^-exp(teta.st) - 
            1)^delta)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * 
            (((1 - p2)^-exp(teta.st) - 1)^(delta - 1) * (delta * 
                ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)))) + 
                ((1 - p1)^-exp(teta.st) - 1)^(delta - 1) * (delta * 
                  ((1 - p1)^-exp(teta.st) * (log((1 - p1)) * 
                    exp(teta.st)))))) * ((1/delta) * (((1 - p2)^-exp(teta.st) - 
            1)^(delta - 1) * (delta * ((1 - p2)^-exp(teta.st) * 
            (log((1 - p2)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 
            1)^(delta - 1) * (delta * ((1 - p1)^-exp(teta.st) * 
            (log((1 - p1)) * exp(teta.st)))))))))




bit1.del2 <-(1 + (((1 - p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^-teta - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((((1 - p1)^-teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * (((1 - p1)^-teta - 
    1)^(exp(delta.st) + 1) * (log(((1 - p1)^-teta - 1)) * exp(delta.st)) + 
    ((1 - p2)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - p2)^-teta - 
        1)) * exp(delta.st)))) - (((1 - p1)^-teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log((((1 - p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - 
    p2)^-teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)))) * ((-1/teta) * ((((1 - p1)^-teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * (((1 - p1)^-teta - 
    1)^(exp(delta.st) + 1) * (log(((1 - p1)^-teta - 1)) * exp(delta.st)) + 
    ((1 - p2)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - p2)^-teta - 
        1)) * exp(delta.st)))) - (((1 - p1)^-teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log((((1 - p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - 
    p2)^-teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)))) + (1 + (((1 - p1)^-teta - 1)^(exp(delta.st) + 1) + 
    ((1 - p2)^-teta - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^((-1/teta) - 1) * ((-1/teta) * (((((1 - p1)^-teta - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1))^(((1/(exp(delta.st) + 1)) - 1) - 1) * (((1/(exp(delta.st) + 
    1)) - 1) * (((1 - p1)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
    p1)^-teta - 1)) * exp(delta.st)) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1) * (log(((1 - p2)^-teta - 1)) * exp(delta.st)))) - (((1 - 
    p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * (log((((1 - p1)^-teta - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * ((1/(exp(delta.st) + 
    1)) * (((1 - p1)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
    p1)^-teta - 1)) * exp(delta.st)) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1) * (log(((1 - p2)^-teta - 1)) * exp(delta.st)))) + (((1 - 
    p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * (((1 - p1)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
    p1)^-teta - 1)) * exp(delta.st)) * (log(((1 - p1)^-teta - 
    1)) * exp(delta.st)) + ((1 - p1)^-teta - 1)^(exp(delta.st) + 
    1) * (log(((1 - p1)^-teta - 1)) * exp(delta.st)) + (((1 - 
    p2)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - p2)^-teta - 
    1)) * exp(delta.st)) * (log(((1 - p2)^-teta - 1)) * exp(delta.st)) + 
    ((1 - p2)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - p2)^-teta - 
        1)) * exp(delta.st)))) - exp(delta.st)/(exp(delta.st) + 
    1)^2 * (((1 - p1)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
    p1)^-teta - 1)) * exp(delta.st)) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1) * (log(((1 - p2)^-teta - 1)) * exp(delta.st)))) - (((((1 - 
    p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * (((1 - p1)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - 
    p1)^-teta - 1)) * exp(delta.st)) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1) * (log(((1 - p2)^-teta - 1)) * exp(delta.st)))) - (((1 - 
    p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)) * (log((((1 - p1)^-teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * (log((((1 - p1)^-teta - 1)^(exp(delta.st) + 1) + 
    ((1 - p2)^-teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)) + (((1 - p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - 
    p2)^-teta - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
    ((((1 - p1)^-teta - 1)^(exp(delta.st) + 1) * (log(((1 - p1)^-teta - 
        1)) * exp(delta.st)) + ((1 - p2)^-teta - 1)^(exp(delta.st) + 
        1) * (log(((1 - p2)^-teta - 1)) * exp(delta.st)))/(((1 - 
        p1)^-teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^-teta - 
        1)^(exp(delta.st) + 1)) * (exp(delta.st)/(exp(delta.st) + 
        1)^2) + log((((1 - p1)^-teta - 1)^(exp(delta.st) + 1) + 
        ((1 - p2)^-teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2 - exp(delta.st) * (2 * (exp(delta.st) * (exp(delta.st) + 
        1)))/((exp(delta.st) + 1)^2)^2)))))




c.copula2.be1del <- ((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p1)^-teta - 1)^delta * log(((1 - p1)^-teta - 1)) + 
        ((1 - p2)^-teta - 1)^delta * log(((1 - p2)^-teta - 1)))) - 
    (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta) * 
        (log((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
            1)^delta)) * (1/delta^2)))) * ((-1/teta) * ((((1 - 
    p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * (((1 - p1)^-teta - 1)^(delta - 1) * (delta * 
    ((1 - p1)^-(teta + 1) * teta)))))) + (1 + (((1 - p1)^-teta - 
    1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * (((1 - 
    p1)^-teta - 1)^delta * log(((1 - p1)^-teta - 1)) + ((1 - 
    p2)^-teta - 1)^delta * log(((1 - p2)^-teta - 1)))) - (((1 - 
    p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
    1) * (log((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)) * (1/delta^2))) * ((1/delta) * (((1 - p1)^-teta - 
    1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 1) * teta)))) + 
    (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * (((1 - p1)^-teta - 1)^(delta - 1) * 
        log(((1 - p1)^-teta - 1)) * (delta * ((1 - p1)^-(teta + 
        1) * teta)) + ((1 - p1)^-teta - 1)^(delta - 1) * ((1 - 
        p1)^-(teta + 1) * teta)) - 1/delta^2 * (((1 - p1)^-teta - 
        1)^(delta - 1) * (delta * ((1 - p1)^-(teta + 1) * teta)))))))*exp(delta.st)


        



        
c.copula2.be2del <- ((1 + (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((((1 - p1)^-teta - 1)^delta + 
    ((1 - p2)^-teta - 1)^delta)^((1/delta) - 1) * ((1/delta) * 
    (((1 - p1)^-teta - 1)^delta * log(((1 - p1)^-teta - 1)) + 
        ((1 - p2)^-teta - 1)^delta * log(((1 - p2)^-teta - 1)))) - 
    (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta) * 
        (log((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
            1)^delta)) * (1/delta^2)))) * ((-1/teta) * ((((1 - 
    p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
    1) * ((1/delta) * (((1 - p2)^-teta - 1)^(delta - 1) * (delta * 
    ((1 - p2)^-(teta + 1) * teta)))))) + (1 + (((1 - p1)^-teta - 
    1)^delta + ((1 - p2)^-teta - 1)^delta)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * (((1 - 
    p1)^-teta - 1)^delta * log(((1 - p1)^-teta - 1)) + ((1 - 
    p2)^-teta - 1)^delta * log(((1 - p2)^-teta - 1)))) - (((1 - 
    p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
    1) * (log((((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 
    1)^delta)) * (1/delta^2))) * ((1/delta) * (((1 - p2)^-teta - 
    1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 1) * teta)))) + 
    (((1 - p1)^-teta - 1)^delta + ((1 - p2)^-teta - 1)^delta)^((1/delta) - 
        1) * ((1/delta) * (((1 - p2)^-teta - 1)^(delta - 1) * 
        log(((1 - p2)^-teta - 1)) * (delta * ((1 - p2)^-(teta + 
        1) * teta)) + ((1 - p2)^-teta - 1)^(delta - 1) * ((1 - 
        p2)^-(teta + 1) * teta)) - 1/delta^2 * (((1 - p2)^-teta - 
        1)^(delta - 1) * (delta * ((1 - p2)^-(teta + 1) * teta)))))))*exp(delta.st)



bit1.thdel <-(1 + (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - 
    p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^((-1/exp(teta.st)) - 1) * ((-1/exp(teta.st)) * ((((1 - 
    p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) * 
    (log(((1 - p1)^-exp(teta.st) - 1)) * exp(delta.st)) + ((1 - 
    p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1) * (log(((1 - p2)^-exp(teta.st) - 
    1)) * exp(delta.st)))) - (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * 
    (exp(delta.st)/(exp(delta.st) + 1)^2)))) * (log((1 + (((1 - 
    p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))) * (exp(teta.st)/exp(teta.st)^2)) + 
    (1 + (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + 
        ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)))^(-1/exp(teta.st)) * (((((1 - p1)^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) * (log(((1 - p1)^-exp(teta.st) - 1)) * exp(delta.st)) + 
        ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1) * (log(((1 - 
            p2)^-exp(teta.st) - 1)) * exp(delta.st)))) - (((1 - 
        p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
        1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log((((1 - 
        p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
        1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2)))/(1 + (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1))) * (exp(teta.st)/exp(teta.st)^2)) - ((1 + (((1 - 
    p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^(((-1/exp(teta.st)) - 
    1) - 1) * (((-1/exp(teta.st)) - 1) * ((((1 - p1)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) * 
    (log(((1 - p1)^-exp(teta.st) - 1)) * exp(delta.st)) + ((1 - 
    p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1) * (log(((1 - p2)^-exp(teta.st) - 
    1)) * exp(delta.st)))) - (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * 
    (exp(delta.st)/(exp(delta.st) + 1)^2)))) * ((-1/exp(teta.st)) * 
    ((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - 
        p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * (((1 - p2)^-exp(teta.st) - 
        1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * 
        ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)))) + 
        ((1 - p1)^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 1) * 
            ((exp(delta.st) + 1) * ((1 - p1)^-exp(teta.st) * 
                (log((1 - p1)) * exp(teta.st)))))))) + (1 + (((1 - 
    p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^((-1/exp(teta.st)) - 
    1) * ((-1/exp(teta.st)) * (((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 
    1)) - 1) - 1) * (((1/(exp(delta.st) + 1)) - 1) * (((1 - p1)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) * (log(((1 - p1)^-exp(teta.st) - 1)) * 
    exp(delta.st)) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) * (log(((1 - p2)^-exp(teta.st) - 1)) * exp(delta.st)))) - 
    (((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - 
        p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * (log((((1 - p1)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * 
        (exp(delta.st)/(exp(delta.st) + 1)^2))) * ((1/(exp(delta.st) + 
    1)) * (((1 - p2)^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 
    1) * ((exp(delta.st) + 1) * ((1 - p2)^-exp(teta.st) * (log((1 - 
    p2)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 1)^((exp(delta.st) + 
    1) - 1) * ((exp(delta.st) + 1) * ((1 - p1)^-exp(teta.st) * 
    (log((1 - p1)) * exp(teta.st)))))) + (((1 - p1)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * (((1 - p2)^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 
    1) * (log(((1 - p2)^-exp(teta.st) - 1)) * exp(delta.st)) * 
    ((exp(delta.st) + 1) * ((1 - p2)^-exp(teta.st) * (log((1 - 
        p2)) * exp(teta.st)))) + ((1 - p2)^-exp(teta.st) - 1)^((exp(delta.st) + 
    1) - 1) * (exp(delta.st) * ((1 - p2)^-exp(teta.st) * (log((1 - 
    p2)) * exp(teta.st)))) + (((1 - p1)^-exp(teta.st) - 1)^((exp(delta.st) + 
    1) - 1) * (log(((1 - p1)^-exp(teta.st) - 1)) * exp(delta.st)) * 
    ((exp(delta.st) + 1) * ((1 - p1)^-exp(teta.st) * (log((1 - 
        p1)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 1)^((exp(delta.st) + 
    1) - 1) * (exp(delta.st) * ((1 - p1)^-exp(teta.st) * (log((1 - 
    p1)) * exp(teta.st)))))) - exp(delta.st)/(exp(delta.st) + 
    1)^2 * (((1 - p2)^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 
    1) * ((exp(delta.st) + 1) * ((1 - p2)^-exp(teta.st) * (log((1 - 
    p2)) * exp(teta.st)))) + ((1 - p1)^-exp(teta.st) - 1)^((exp(delta.st) + 
    1) - 1) * ((exp(delta.st) + 1) * ((1 - p1)^-exp(teta.st) * 
    (log((1 - p1)) * exp(teta.st)))))))))


}


if(BivD=="BB1.270"){
  c.copula.be1 <- 1 - (1 + ((p1^(teta) - 1)^(-delta) + ((1 - p2)^(teta) - 1)^(-delta))^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((p1^(teta) - 1)^(-delta) + ((1 - p2)^(teta) - 
    1)^(-delta))^((-1/delta) - 1) * ((-1/delta) * ((p1^(teta) - 
    1)^((-delta) - 1) * ((-delta) * (p1^((teta) - 1) * (teta)))))))



 
  c.copula.be2 <- -((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - p2)^teta - 
    1)^-(delta + 1) * (delta * ((1 - p2)^(teta - 1) * teta)))))))
    
    
  c.copula.theta <- ((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^(1/teta) * 
    (log((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))) * 
        (1/teta^2)) + (1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^(-1/delta))^((1/teta) - 1) * ((1/teta) * (((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 1) * 
    ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 1) * (delta * 
        ((1 - p2)^teta * log((1 - p2)))) + (p1^teta - 1)^-(delta + 
        1) * (delta * (p1^teta * log(p1))))))))*(-exp(teta.st))


  c.copula.delta <- (-((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^(-1/delta) * (log(((p1^teta - 1)^-delta + ((1 - 
    p2)^teta - 1)^-delta)) * (1/delta^2)) - ((p1^teta - 1)^-delta + 
    ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * 
    (((1 - p2)^teta - 1)^-delta * log(((1 - p2)^teta - 1)) + 
        (p1^teta - 1)^-delta * log((p1^teta - 1))))))))*(-exp(delta.st))



 c.copula2.be1 <-  (1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p1^teta - 1)^-(delta + 
    1) * (delta * (p1^((teta - 1) - 1) * (teta - 1) * teta)) - 
    (p1^teta - 1)^-(delta + 1 + 1) * ((delta + 1) * (p1^(teta - 
        1) * teta)) * (delta * (p1^(teta - 1) * teta)))) - ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^(((-1/delta) - 1) - 
    1) * (((-1/delta) - 1) * ((p1^teta - 1)^-(delta + 1) * (delta * 
    (p1^(teta - 1) * teta)))) * ((-1/delta) * ((p1^teta - 1)^-(delta + 
    1) * (delta * (p1^(teta - 1) * teta)))))) - (1 + ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (((p1^teta - 1)^-delta + ((1 - 
    p2)^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p1^teta - 
    1)^-(delta + 1) * (delta * (p1^(teta - 1) * teta)))))) * 
    ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((p1^teta - 1)^-(delta + 1) * (delta * 
        (p1^(teta - 1) * teta))))))


                  
 c.copula2.be2 <- -((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (((p1^teta - 1)^-delta + ((1 - 
    p2)^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - 
    p2)^teta - 1)^-(delta + 1) * (delta * ((1 - p2)^(teta - 1) * 
    teta)))))) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - p2)^teta - 
    1)^-(delta + 1) * (delta * ((1 - p2)^(teta - 1) * teta)))))) + 
    (1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
        1) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
        1)^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * 
        (((1 - p2)^teta - 1)^-(delta + 1) * (delta * ((1 - p2)^(teta - 
            1) * teta)))) * ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 
        1) * (delta * ((1 - p2)^(teta - 1) * teta)))) + ((p1^teta - 
        1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 1 + 
        1) * ((delta + 1) * ((1 - p2)^(teta - 1) * teta)) * (delta * 
        ((1 - p2)^(teta - 1) * teta)) - ((1 - p2)^teta - 1)^-(delta + 
        1) * (delta * ((1 - p2)^((teta - 1) - 1) * (teta - 1) * 
        teta)))))))

                  
                 



c.copula2.be1be2 <- (1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (((p1^teta - 1)^-delta + ((1 - 
    p2)^teta - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - 
    p2)^teta - 1)^-(delta + 1) * (delta * ((1 - p2)^(teta - 1) * 
    teta)))))) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p1^teta - 1)^-(delta + 
    1) * (delta * (p1^(teta - 1) * teta)))))) + (1 + ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * (((1 - 
    p2)^teta - 1)^-(delta + 1) * (delta * ((1 - p2)^(teta - 1) * 
    teta)))) * ((-1/delta) * ((p1^teta - 1)^-(delta + 1) * (delta * 
    (p1^(teta - 1) * teta))))))



c.copula2.be1th <-((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p1^teta - 1)^-(delta + 
    1) * (delta * (p1^(teta - 1) * log(p1) * teta + p1^(teta - 
    1))) - (p1^teta - 1)^-(delta + 1 + 1) * ((delta + 1) * (p1^teta * 
    log(p1))) * (delta * (p1^(teta - 1) * teta)))) - ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^(((-1/delta) - 1) - 
    1) * (((-1/delta) - 1) * (((1 - p2)^teta - 1)^-(delta + 1) * 
    (delta * ((1 - p2)^teta * log((1 - p2)))) + (p1^teta - 1)^-(delta + 
    1) * (delta * (p1^teta * log(p1))))) * ((-1/delta) * ((p1^teta - 
    1)^-(delta + 1) * (delta * (p1^(teta - 1) * teta))))) - 1/teta^2 * 
    (((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((p1^teta - 1)^-(delta + 1) * (delta * 
        (p1^(teta - 1) * teta)))))) - ((1 + ((p1^teta - 1)^-delta + 
    ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 1) * 
    (log((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))) * 
        (1/teta^2)) + (1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 
    1) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 1) * (delta * 
    ((1 - p2)^teta * log((1 - p2)))) + (p1^teta - 1)^-(delta + 
    1) * (delta * (p1^teta * log(p1)))))))) * ((1/teta) * (((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 1) * 
    ((-1/delta) * ((p1^teta - 1)^-(delta + 1) * (delta * (p1^(teta - 
        1) * teta)))))))*(-exp(teta.st))


        

c.copula2.be2th <-(-((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * (((1 - p2)^teta - 
    1)^-(delta + 1) * (delta * ((1 - p2)^(teta - 1) * log((1 - 
    p2)) * teta + (1 - p2)^(teta - 1))) - ((1 - p2)^teta - 1)^-(delta + 
    1 + 1) * ((delta + 1) * ((1 - p2)^teta * log((1 - p2)))) * 
    (delta * ((1 - p2)^(teta - 1) * teta)))) - ((p1^teta - 1)^-delta + 
    ((1 - p2)^teta - 1)^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
    1) * (((1 - p2)^teta - 1)^-(delta + 1) * (delta * ((1 - p2)^teta * 
    log((1 - p2)))) + (p1^teta - 1)^-(delta + 1) * (delta * (p1^teta * 
    log(p1))))) * ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 
    1) * (delta * ((1 - p2)^(teta - 1) * teta))))) - 1/teta^2 * 
    (((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 1) * 
        (delta * ((1 - p2)^(teta - 1) * teta)))))) - ((1 + ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * (log((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))) * 
    (1/teta^2)) + (1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 
    1) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 1) * (delta * 
    ((1 - p2)^teta * log((1 - p2)))) + (p1^teta - 1)^-(delta + 
    1) * (delta * (p1^teta * log(p1)))))))) * ((1/teta) * (((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 1) * 
    ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 1) * (delta * 
        ((1 - p2)^(teta - 1) * teta))))))))*(-exp(teta.st))






bit1.th2 <--(((1 + ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
    1)^-delta)^(-1/delta))^(((-1/exp(teta.st)) - 1) - 1) * (((-1/exp(teta.st)) - 
    1) * (((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
    1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p1^-exp(teta.st) - 
    1)^-(delta + 1) * (delta * (p1^-exp(teta.st) * (log(p1) * 
    exp(teta.st)))) + ((1 - p2)^-exp(teta.st) - 1)^-(delta + 
    1) * (delta * ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * 
    exp(teta.st)))))))) + (1 + ((p1^-exp(teta.st) - 1)^-delta + 
    ((1 - p2)^-exp(teta.st) - 1)^-delta)^(-1/delta))^((-1/exp(teta.st)) - 
    1) * (log((1 + ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
    1)^-delta)^(-1/delta))) * (exp(teta.st)/exp(teta.st)^2))) * 
    ((-1/exp(teta.st)) * (((p1^-exp(teta.st) - 1)^-delta + ((1 - 
        p2)^-exp(teta.st) - 1)^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((p1^-exp(teta.st) - 1)^-(delta + 1) * (delta * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st)))) + ((1 - p2)^-exp(teta.st) - 
            1)^-(delta + 1) * (delta * ((1 - p2)^-exp(teta.st) * 
            (log((1 - p2)) * exp(teta.st)))))))) + (1 + ((p1^-exp(teta.st) - 
    1)^-delta + ((1 - p2)^-exp(teta.st) - 1)^-delta)^(-1/delta))^((-1/exp(teta.st)) - 
    1) * (exp(teta.st)/exp(teta.st)^2 * (((p1^-exp(teta.st) - 
    1)^-delta + ((1 - p2)^-exp(teta.st) - 1)^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((p1^-exp(teta.st) - 1)^-(delta + 1) * 
    (delta * (p1^-exp(teta.st) * (log(p1) * exp(teta.st)))) + 
    ((1 - p2)^-exp(teta.st) - 1)^-(delta + 1) * (delta * ((1 - 
        p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st))))))) + 
    (-1/exp(teta.st)) * (((p1^-exp(teta.st) - 1)^-delta + ((1 - 
        p2)^-exp(teta.st) - 1)^-delta)^(((-1/delta) - 1) - 1) * 
        (((-1/delta) - 1) * ((p1^-exp(teta.st) - 1)^-(delta + 
            1) * (delta * (p1^-exp(teta.st) * (log(p1) * exp(teta.st)))) + 
            ((1 - p2)^-exp(teta.st) - 1)^-(delta + 1) * (delta * 
                ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)))))) * 
        ((-1/delta) * ((p1^-exp(teta.st) - 1)^-(delta + 1) * 
            (delta * (p1^-exp(teta.st) * (log(p1) * exp(teta.st)))) + 
            ((1 - p2)^-exp(teta.st) - 1)^-(delta + 1) * (delta * 
                ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)))))) + 
        ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
            1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p1^-exp(teta.st) - 
            1)^-(delta + 1 + 1) * ((delta + 1) * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st)))) * (delta * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st)))) + (p1^-exp(teta.st) - 
            1)^-(delta + 1) * (delta * (p1^-exp(teta.st) * (log(p1) * 
            exp(teta.st)) - p1^-exp(teta.st) * (log(p1) * exp(teta.st)) * 
            (log(p1) * exp(teta.st)))) + (((1 - p2)^-exp(teta.st) - 
            1)^-(delta + 1 + 1) * ((delta + 1) * ((1 - p2)^-exp(teta.st) * 
            (log((1 - p2)) * exp(teta.st)))) * (delta * ((1 - 
            p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)))) + 
            ((1 - p2)^-exp(teta.st) - 1)^-(delta + 1) * (delta * 
                ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)) - 
                  (1 - p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)) * 
                    (log((1 - p2)) * exp(teta.st))))))))) + (((1 + 
    ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
        1)^-delta)^(-1/delta))^((-1/exp(teta.st)) - 1) * ((-1/exp(teta.st)) * 
    (((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
        1)^-delta)^((-1/delta) - 1) * ((-1/delta) * ((p1^-exp(teta.st) - 
        1)^-(delta + 1) * (delta * (p1^-exp(teta.st) * (log(p1) * 
        exp(teta.st)))) + ((1 - p2)^-exp(teta.st) - 1)^-(delta + 
        1) * (delta * ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * 
        exp(teta.st)))))))) + (1 + ((p1^-exp(teta.st) - 1)^-delta + 
    ((1 - p2)^-exp(teta.st) - 1)^-delta)^(-1/delta))^(-1/exp(teta.st)) * 
    (log((1 + ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
        1)^-delta)^(-1/delta))) * (exp(teta.st)/exp(teta.st)^2))) * 
    (log((1 + ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
        1)^-delta)^(-1/delta))) * (exp(teta.st)/exp(teta.st)^2)) + 
    (1 + ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
        1)^-delta)^(-1/delta))^(-1/exp(teta.st)) * (((p1^-exp(teta.st) - 
        1)^-delta + ((1 - p2)^-exp(teta.st) - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((p1^-exp(teta.st) - 1)^-(delta + 
        1) * (delta * (p1^-exp(teta.st) * (log(p1) * exp(teta.st)))) + 
        ((1 - p2)^-exp(teta.st) - 1)^-(delta + 1) * (delta * 
            ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st))))))/(1 + 
        ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
            1)^-delta)^(-1/delta)) * (exp(teta.st)/exp(teta.st)^2) + 
        log((1 + ((p1^-exp(teta.st) - 1)^-delta + ((1 - p2)^-exp(teta.st) - 
            1)^-delta)^(-1/delta))) * (exp(teta.st)/exp(teta.st)^2 - 
            exp(teta.st) * (2 * (exp(teta.st) * exp(teta.st)))/(exp(teta.st)^2)^2))))






bit1.del2 <--((1 + ((p1^teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^teta - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (((p1^teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p1^teta - 1)^(exp(delta.st) + 
    1) * (log((p1^teta - 1)) * exp(delta.st)) + ((1 - p2)^teta - 
    1)^(exp(delta.st) + 1) * (log(((1 - p2)^teta - 1)) * exp(delta.st)))) - 
    ((p1^teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log(((p1^teta - 1)^(exp(delta.st) + 
        1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2)))) * ((1/teta) * (((p1^teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p1^teta - 1)^(exp(delta.st) + 
    1) * (log((p1^teta - 1)) * exp(delta.st)) + ((1 - p2)^teta - 
    1)^(exp(delta.st) + 1) * (log(((1 - p2)^teta - 1)) * exp(delta.st)))) - 
    ((p1^teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log(((p1^teta - 1)^(exp(delta.st) + 
        1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2)))) + (1 + ((p1^teta - 1)^(exp(delta.st) + 1) + 
    ((1 - p2)^teta - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^((1/teta) - 1) * ((1/teta) * ((((p1^teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 
    1)) - 1) - 1) * (((1/(exp(delta.st) + 1)) - 1) * ((p1^teta - 
    1)^(exp(delta.st) + 1) * (log((p1^teta - 1)) * exp(delta.st)) + 
    ((1 - p2)^teta - 1)^(exp(delta.st) + 1) * (log(((1 - p2)^teta - 
        1)) * exp(delta.st)))) - ((p1^teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * (log(((p1^teta - 1)^(exp(delta.st) + 1) + ((1 - 
    p2)^teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * ((1/(exp(delta.st) + 1)) * ((p1^teta - 1)^(exp(delta.st) + 
    1) * (log((p1^teta - 1)) * exp(delta.st)) + ((1 - p2)^teta - 
    1)^(exp(delta.st) + 1) * (log(((1 - p2)^teta - 1)) * exp(delta.st)))) + 
    ((p1^teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * ((p1^teta - 1)^(exp(delta.st) + 1) * (log((p1^teta - 
        1)) * exp(delta.st)) * (log((p1^teta - 1)) * exp(delta.st)) + 
        (p1^teta - 1)^(exp(delta.st) + 1) * (log((p1^teta - 1)) * 
            exp(delta.st)) + (((1 - p2)^teta - 1)^(exp(delta.st) + 
        1) * (log(((1 - p2)^teta - 1)) * exp(delta.st)) * (log(((1 - 
        p2)^teta - 1)) * exp(delta.st)) + ((1 - p2)^teta - 1)^(exp(delta.st) + 
        1) * (log(((1 - p2)^teta - 1)) * exp(delta.st)))) - exp(delta.st)/(exp(delta.st) + 
        1)^2 * ((p1^teta - 1)^(exp(delta.st) + 1) * (log((p1^teta - 
        1)) * exp(delta.st)) + ((1 - p2)^teta - 1)^(exp(delta.st) + 
        1) * (log(((1 - p2)^teta - 1)) * exp(delta.st)))) - ((((p1^teta - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((p1^teta - 1)^(exp(delta.st) + 1) * (log((p1^teta - 
    1)) * exp(delta.st)) + ((1 - p2)^teta - 1)^(exp(delta.st) + 
    1) * (log(((1 - p2)^teta - 1)) * exp(delta.st)))) - ((p1^teta - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)) * (log(((p1^teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * (log(((p1^teta - 1)^(exp(delta.st) + 1) + ((1 - 
    p2)^teta - 1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)) + ((p1^teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^teta - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (((p1^teta - 
    1)^(exp(delta.st) + 1) * (log((p1^teta - 1)) * exp(delta.st)) + 
    ((1 - p2)^teta - 1)^(exp(delta.st) + 1) * (log(((1 - p2)^teta - 
        1)) * exp(delta.st)))/((p1^teta - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^teta - 1)^(exp(delta.st) + 1)) * (exp(delta.st)/(exp(delta.st) + 
    1)^2) + log(((p1^teta - 1)^(exp(delta.st) + 1) + ((1 - p2)^teta - 
    1)^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2 - exp(delta.st) * (2 * (exp(delta.st) * (exp(delta.st) + 
    1)))/((exp(delta.st) + 1)^2)^2))))))




c.copula2.be1del <- ((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (((p1^teta - 1)^-delta + ((1 - 
    p2)^teta - 1)^-delta)^(-1/delta) * (log(((p1^teta - 1)^-delta + 
    ((1 - p2)^teta - 1)^-delta)) * (1/delta^2)) - ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 1) * 
    ((-1/delta) * (((1 - p2)^teta - 1)^-delta * log(((1 - p2)^teta - 
        1)) + (p1^teta - 1)^-delta * log((p1^teta - 1)))))) * 
    ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((p1^teta - 1)^-(delta + 1) * (delta * 
        (p1^(teta - 1) * teta)))))) + (1 + ((p1^teta - 1)^-delta + 
    ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 1) * 
    ((1/teta) * ((((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
        1) * (log(((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)) * 
        (1/delta^2)) - ((p1^teta - 1)^-delta + ((1 - p2)^teta - 
        1)^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * 
        (((1 - p2)^teta - 1)^-delta * log(((1 - p2)^teta - 1)) + 
            (p1^teta - 1)^-delta * log((p1^teta - 1))))) * ((-1/delta) * 
        ((p1^teta - 1)^-(delta + 1) * (delta * (p1^(teta - 1) * 
            teta)))) + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 
        1)^-delta)^((-1/delta) - 1) * (1/delta^2 * ((p1^teta - 
        1)^-(delta + 1) * (delta * (p1^(teta - 1) * teta))) + 
        (-1/delta) * ((p1^teta - 1)^-(delta + 1) * (p1^(teta - 
            1) * teta) - (p1^teta - 1)^-(delta + 1) * log((p1^teta - 
            1)) * (delta * (p1^(teta - 1) * teta)))))))*(-exp(delta.st))


        
        


        
c.copula2.be2del <- (-((1 + ((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (((p1^teta - 1)^-delta + ((1 - 
    p2)^teta - 1)^-delta)^(-1/delta) * (log(((p1^teta - 1)^-delta + 
    ((1 - p2)^teta - 1)^-delta)) * (1/delta^2)) - ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 1) * 
    ((-1/delta) * (((1 - p2)^teta - 1)^-delta * log(((1 - p2)^teta - 
        1)) + (p1^teta - 1)^-delta * log((p1^teta - 1)))))) * 
    ((1/teta) * (((p1^teta - 1)^-delta + ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 
        1) * ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 1) * 
        (delta * ((1 - p2)^(teta - 1) * teta)))))) + (1 + ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((p1^teta - 1)^-delta + ((1 - p2)^teta - 
    1)^-delta)^((-1/delta) - 1) * (log(((p1^teta - 1)^-delta + 
    ((1 - p2)^teta - 1)^-delta)) * (1/delta^2)) - ((p1^teta - 
    1)^-delta + ((1 - p2)^teta - 1)^-delta)^(((-1/delta) - 1) - 
    1) * (((-1/delta) - 1) * (((1 - p2)^teta - 1)^-delta * log(((1 - 
    p2)^teta - 1)) + (p1^teta - 1)^-delta * log((p1^teta - 1))))) * 
    ((-1/delta) * (((1 - p2)^teta - 1)^-(delta + 1) * (delta * 
        ((1 - p2)^(teta - 1) * teta)))) + ((p1^teta - 1)^-delta + 
    ((1 - p2)^teta - 1)^-delta)^((-1/delta) - 1) * (1/delta^2 * 
    (((1 - p2)^teta - 1)^-(delta + 1) * (delta * ((1 - p2)^(teta - 
        1) * teta))) + (-1/delta) * (((1 - p2)^teta - 1)^-(delta + 
    1) * ((1 - p2)^(teta - 1) * teta) - ((1 - p2)^teta - 1)^-(delta + 
    1) * log(((1 - p2)^teta - 1)) * (delta * ((1 - p2)^(teta - 
    1) * teta)))))))
)*(-exp(delta.st))




bit1.thdel <--((1 + ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
    1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^((-1/exp(teta.st)) - 
    1) * ((-1/exp(teta.st)) * (((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) * (log((p1^-exp(teta.st) - 1)) * exp(delta.st)) + 
    ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1) * (log(((1 - 
        p2)^-exp(teta.st) - 1)) * exp(delta.st)))) - ((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)) * (log(((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
    1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))) * 
    (exp(delta.st)/(exp(delta.st) + 1)^2)))) * (log((1 + ((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))) * (exp(teta.st)/exp(teta.st)^2)) + 
    (1 + ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - 
        p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)))^(-1/exp(teta.st)) * ((((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p1^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) * (log((p1^-exp(teta.st) - 1)) * 
        exp(delta.st)) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) * (log(((1 - p2)^-exp(teta.st) - 1)) * exp(delta.st)))) - 
        ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
            1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
            (log(((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + 
                ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
                  1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)))/(1 + 
        ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
            1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        (exp(teta.st)/exp(teta.st)^2)) - ((1 + ((p1^-exp(teta.st) - 
    1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^(((-1/exp(teta.st)) - 1) - 1) * 
    (((-1/exp(teta.st)) - 1) * (((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((p1^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) * (log((p1^-exp(teta.st) - 1)) * 
        exp(delta.st)) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) * (log(((1 - p2)^-exp(teta.st) - 1)) * exp(delta.st)))) - 
        ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 
            1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
            (log(((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + 
                ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
                  1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)))) * 
    ((-1/exp(teta.st)) * (((p1^-exp(teta.st) - 1)^(exp(delta.st) + 
        1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * (((1 - p2)^-exp(teta.st) - 
        1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * 
        ((1 - p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)))) + 
        (p1^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
            1) * (p1^-exp(teta.st) * (log(p1) * exp(teta.st)))))))) + 
    (1 + ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) + ((1 - 
        p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)))^((-1/exp(teta.st)) - 1) * ((-1/exp(teta.st)) * ((((p1^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1))^(((1/(exp(delta.st) + 1)) - 1) - 1) * (((1/(exp(delta.st) + 
        1)) - 1) * ((p1^-exp(teta.st) - 1)^(exp(delta.st) + 1) * 
        (log((p1^-exp(teta.st) - 1)) * exp(delta.st)) + ((1 - 
        p2)^-exp(teta.st) - 1)^(exp(delta.st) + 1) * (log(((1 - 
        p2)^-exp(teta.st) - 1)) * exp(delta.st)))) - ((p1^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * (log(((p1^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * ((1/(exp(delta.st) + 
        1)) * (((1 - p2)^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * ((exp(delta.st) + 1) * ((1 - p2)^-exp(teta.st) * 
        (log((1 - p2)) * exp(teta.st)))) + (p1^-exp(teta.st) - 
        1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * 
        (p1^-exp(teta.st) * (log(p1) * exp(teta.st)))))) + ((p1^-exp(teta.st) - 
        1)^(exp(delta.st) + 1) + ((1 - p2)^-exp(teta.st) - 1)^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * (((1 - p2)^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * (log(((1 - p2)^-exp(teta.st) - 1)) * exp(delta.st)) * 
        ((exp(delta.st) + 1) * ((1 - p2)^-exp(teta.st) * (log((1 - 
            p2)) * exp(teta.st)))) + ((1 - p2)^-exp(teta.st) - 
        1)^((exp(delta.st) + 1) - 1) * (exp(delta.st) * ((1 - 
        p2)^-exp(teta.st) * (log((1 - p2)) * exp(teta.st)))) + 
        ((p1^-exp(teta.st) - 1)^((exp(delta.st) + 1) - 1) * (log((p1^-exp(teta.st) - 
            1)) * exp(delta.st)) * ((exp(delta.st) + 1) * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st)))) + (p1^-exp(teta.st) - 
            1)^((exp(delta.st) + 1) - 1) * (exp(delta.st) * (p1^-exp(teta.st) * 
            (log(p1) * exp(teta.st)))))) - exp(delta.st)/(exp(delta.st) + 
        1)^2 * (((1 - p2)^-exp(teta.st) - 1)^((exp(delta.st) + 
        1) - 1) * ((exp(delta.st) + 1) * ((1 - p2)^-exp(teta.st) * 
        (log((1 - p2)) * exp(teta.st)))) + (p1^-exp(teta.st) - 
        1)^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * 
        (p1^-exp(teta.st) * (log(p1) * exp(teta.st))))))))))


}



if(BivD=="BB6.0"){
   
  c.copula.be1 <- exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 1) * 
    (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta))))))))

 
  c.copula.be2 <- exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 1) * 
    (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta))))))))

  c.copula.theta <- (-(exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^(1/teta) * 
    (log(exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - 
        p2)^teta))^delta)^(1/delta))) * (1/teta^2)) + exp(-((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 1) * 
    (delta * ((1 - p1)^teta * log((1 - p1))/(1 - (1 - p1)^teta))) + 
    (-log(1 - (1 - p2)^teta))^(delta - 1) * (delta * ((1 - p2)^teta * 
        log((1 - p2))/(1 - (1 - p2)^teta))))))))))*exp(teta.st)



  c.copula.delta <- (-(exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^delta * log((-log(1 - 
    (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^delta * log((-log(1 - 
    (1 - p2)^teta))))) - ((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
    (1 - p2)^teta))^delta)^(1/delta) * (log(((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)) * (1/delta^2)))))))*exp(delta.st)




 c.copula2.be1 <-  exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 1) * 
    (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta)))))))) * 
    ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
        (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - (1 - 
        p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
        1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta)))))))) + 
    exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
        1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
        1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta)))))) * 
        (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
            1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - 
            p1)^teta)))))) - exp(-((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
        1) * (delta * ((1 - p1)^((teta - 1) - 1) * (teta - 1) * 
        teta/(1 - (1 - p1)^teta) + (1 - p1)^(teta - 1) * teta * 
        ((1 - p1)^(teta - 1) * teta)/(1 - (1 - p1)^teta)^2)) + 
        (-log(1 - (1 - p1)^teta))^((delta - 1) - 1) * ((delta - 
            1) * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta))) * 
            (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta))))) + 
        ((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(((1/delta) - 
            1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - p1)^teta))^(delta - 
            1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - 
            p1)^teta))))) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
            1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - 
            p1)^teta))))))))


                  
 c.copula2.be2 <- exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 1) * 
    (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta)))))))) * 
    ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
        (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - (1 - 
        p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
        1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta)))))))) + 
    exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
        1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
        1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta)))))) * 
        (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
            1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - 
            p2)^teta)))))) - exp(-((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
        1) * (delta * ((1 - p2)^((teta - 1) - 1) * (teta - 1) * 
        teta/(1 - (1 - p2)^teta) + (1 - p2)^(teta - 1) * teta * 
        ((1 - p2)^(teta - 1) * teta)/(1 - (1 - p2)^teta)^2)) + 
        (-log(1 - (1 - p2)^teta))^((delta - 1) - 1) * ((delta - 
            1) * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta))) * 
            (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta))))) + 
        ((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(((1/delta) - 
            1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - p2)^teta))^(delta - 
            1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - 
            p2)^teta))))) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
            1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - 
            p2)^teta))))))))


                  
                  
             


c.copula2.be1be2 <- exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 1) * 
    (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta)))))))) * 
    ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
        (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - (1 - 
        p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
        1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta)))))))) + 
    exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
        1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
        1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta)))))) * 
        (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
            1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - 
            p1)^teta)))))) - exp(-((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - p2)^teta))^(delta - 
        1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta))))) * 
        ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 1) * 
            (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta))))))))




c.copula2.be1th <-(exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - p1)^teta))^(delta - 
    1) * (delta * ((1 - p1)^teta * log((1 - p1))/(1 - (1 - p1)^teta))) + 
    (-log(1 - (1 - p2)^teta))^(delta - 1) * (delta * ((1 - p2)^teta * 
        log((1 - p2))/(1 - (1 - p2)^teta))))) * ((1/delta) * 
    ((-log(1 - (1 - p1)^teta))^(delta - 1) * (delta * ((1 - p1)^(teta - 
        1) * teta/(1 - (1 - p1)^teta))))) + ((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
    ((-log(1 - (1 - p1)^teta))^((delta - 1) - 1) * ((delta - 
        1) * ((1 - p1)^teta * log((1 - p1))/(1 - (1 - p1)^teta))) * 
        (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta))) + 
        (-log(1 - (1 - p1)^teta))^(delta - 1) * (delta * (((1 - 
            p1)^(teta - 1) * log((1 - p1)) * teta + (1 - p1)^(teta - 
            1))/(1 - (1 - p1)^teta) + (1 - p1)^(teta - 1) * teta * 
            ((1 - p1)^teta * log((1 - p1)))/(1 - (1 - p1)^teta)^2))))) - 
    exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
            1) * (delta * ((1 - p1)^teta * log((1 - p1))/(1 - 
            (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^(delta - 
            1) * (delta * ((1 - p2)^teta * log((1 - p2))/(1 - 
            (1 - p2)^teta)))))) * (((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - p1)^teta))^(delta - 1) * (delta * ((1 - 
            p1)^(teta - 1) * teta/(1 - (1 - p1)^teta))))))) - 
    1/teta^2 * (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
        (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - (1 - 
        p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
        1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta)))))))) - 
    (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
        1) * (log(exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
        (1 - p2)^teta))^delta)^(1/delta))) * (1/teta^2)) + exp(-((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
        1) * (delta * ((1 - p1)^teta * log((1 - p1))/(1 - (1 - 
        p1)^teta))) + (-log(1 - (1 - p2)^teta))^(delta - 1) * 
        (delta * ((1 - p2)^teta * log((1 - p2))/(1 - (1 - p2)^teta))))))))) * 
        ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
            (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
            (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
            1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - 
            p1)^teta)))))))))*exp(teta.st)


        


c.copula2.be2th <-(exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - p1)^teta))^(delta - 
    1) * (delta * ((1 - p1)^teta * log((1 - p1))/(1 - (1 - p1)^teta))) + 
    (-log(1 - (1 - p2)^teta))^(delta - 1) * (delta * ((1 - p2)^teta * 
        log((1 - p2))/(1 - (1 - p2)^teta))))) * ((1/delta) * 
    ((-log(1 - (1 - p2)^teta))^(delta - 1) * (delta * ((1 - p2)^(teta - 
        1) * teta/(1 - (1 - p2)^teta))))) + ((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
    ((-log(1 - (1 - p2)^teta))^((delta - 1) - 1) * ((delta - 
        1) * ((1 - p2)^teta * log((1 - p2))/(1 - (1 - p2)^teta))) * 
        (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta))) + 
        (-log(1 - (1 - p2)^teta))^(delta - 1) * (delta * (((1 - 
            p2)^(teta - 1) * log((1 - p2)) * teta + (1 - p2)^(teta - 
            1))/(1 - (1 - p2)^teta) + (1 - p2)^(teta - 1) * teta * 
            ((1 - p2)^teta * log((1 - p2)))/(1 - (1 - p2)^teta)^2))))) - 
    exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
            1) * (delta * ((1 - p1)^teta * log((1 - p1))/(1 - 
            (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^(delta - 
            1) * (delta * ((1 - p2)^teta * log((1 - p2))/(1 - 
            (1 - p2)^teta)))))) * (((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - p2)^teta))^(delta - 1) * (delta * ((1 - 
            p2)^(teta - 1) * teta/(1 - (1 - p2)^teta))))))) - 
    1/teta^2 * (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
        (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - (1 - 
        p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
        1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta)))))))) - 
    (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
        1) * (log(exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
        (1 - p2)^teta))^delta)^(1/delta))) * (1/teta^2)) + exp(-((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
        1) * (delta * ((1 - p1)^teta * log((1 - p1))/(1 - (1 - 
        p1)^teta))) + (-log(1 - (1 - p2)^teta))^(delta - 1) * 
        (delta * ((1 - p2)^teta * log((1 - p2))/(1 - (1 - p2)^teta))))))))) * 
        ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
            (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
            (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
            1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - 
            p2)^teta)))))))))*exp(teta.st)





bit1.th2 <--(exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta))^(1/(exp(teta.st) + 
    1)) * (log(exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * (2 * 
        (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
        1)^2)^2) - exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
    (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(delta - 
            1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
            p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(delta - 
            1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
            p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1)))))))/exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) - (exp(-((-log(1 - (1 - 
    p1)^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^delta)^(1/delta))^(1/(exp(teta.st) + 1)) * (log(exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^delta)^(1/delta))) * (exp(teta.st)/(exp(teta.st) + 1)^2)) + 
    exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta))^((1/(exp(teta.st) + 
        1)) - 1) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - 
        (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
        p2)^(exp(teta.st) + 1)))^delta)^(1/delta)) * (((-log(1 - 
        (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
        p2)^(exp(teta.st) + 1)))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(delta - 1) * 
            (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
                p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
                1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(delta - 
            1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
            p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1)))))))))) * (log(exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) + (exp(-((-log(1 - (1 - 
    p1)^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^delta)^(1/delta))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 
    1)) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
    (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^delta)^(((1/delta) - 1) - 
        1) * (((1/delta) - 1) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(delta - 
        1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
        p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 1)))))) * 
        ((1/delta) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(delta - 
            1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
            p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(delta - 
            1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
            p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1)))))) + ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + 
        (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^((delta - 1) - 1) * ((delta - 1) * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
        1)))) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
        p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 1)))) + 
        (-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(delta - 1) * 
            (delta * (((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
                p1)) * exp(teta.st)) * (log((1 - p1)) * exp(teta.st)) + 
                (1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
                  exp(teta.st)))/(1 - (1 - p1)^(exp(teta.st) + 
                1)) + (1 - p1)^(exp(teta.st) + 1) * (log((1 - 
                p1)) * exp(teta.st)) * ((1 - p1)^(exp(teta.st) + 
                1) * (log((1 - p1)) * exp(teta.st)))/(1 - (1 - 
                p1)^(exp(teta.st) + 1))^2)) + ((-log(1 - (1 - 
        p2)^(exp(teta.st) + 1)))^((delta - 1) - 1) * ((delta - 
        1) * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
        exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 1)))) * (delta * 
        ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
            (1 - p2)^(exp(teta.st) + 1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * (((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st)) * (log((1 - p2)) * 
        exp(teta.st)) + (1 - p2)^(exp(teta.st) + 1) * (log((1 - 
        p2)) * exp(teta.st)))/(1 - (1 - p2)^(exp(teta.st) + 1)) + 
        (1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)) * 
            ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)))/(1 - 
            (1 - p2)^(exp(teta.st) + 1))^2)))))) - exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^delta)^(1/delta)) * (((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(delta - 
    1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
    exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 1)))) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(delta - 1) * (delta * ((1 - 
    p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
    (1 - p2)^(exp(teta.st) + 1))))))) * (((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(delta - 
    1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
    exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 1)))) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(delta - 1) * (delta * ((1 - 
    p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
    (1 - p2)^(exp(teta.st) + 1)))))))) - exp(teta.st)/(exp(teta.st) + 
    1)^2 * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
    (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(delta - 
            1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
            p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(delta - 
            1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
            p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1))))))))) - (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * (log(exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) + exp(-((-log(1 - (1 - 
    p1)^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^delta)^(1/delta))^(((1/(exp(teta.st) + 1)) - 1) - 1) * 
    (((1/(exp(teta.st) + 1)) - 1) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
        (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - 
            (1 - p2)^(exp(teta.st) + 1)))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
            1)))^(delta - 1) * (delta * ((1 - p1)^(exp(teta.st) + 
            1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(delta - 
            1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
            p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1)))))))))) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^delta)^(1/delta)) * (((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(delta - 
    1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
    exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 1)))) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(delta - 1) * (delta * ((1 - 
    p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
    (1 - p2)^(exp(teta.st) + 1)))))))))))




bit1.del2 <- -(exp(-((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^((1/teta) - 1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * ((((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 
    1)) - 1) - 1) * (((1/(exp(delta.st) + 1)) - 1) * ((-log(1 - 
    (1 - p1)^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    p1)^teta))) * exp(delta.st)) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p2)^teta))) * exp(delta.st)))) - 
    ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - p2)^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * (log(((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))) * 
        (exp(delta.st)/(exp(delta.st) + 1)^2))) * ((1/(exp(delta.st) + 
    1)) * ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - p1)^teta))) * exp(delta.st)) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p2)^teta))) * exp(delta.st)))) + 
    ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - p2)^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - 
        p1)^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
        p1)^teta))) * exp(delta.st)) * (log((-log(1 - (1 - p1)^teta))) * 
        exp(delta.st)) + (-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
        1) * (log((-log(1 - (1 - p1)^teta))) * exp(delta.st)) + 
        ((-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
            (1 - p2)^teta))) * exp(delta.st)) * (log((-log(1 - 
            (1 - p2)^teta))) * exp(delta.st)) + (-log(1 - (1 - 
            p2)^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
            p2)^teta))) * exp(delta.st)))) - exp(delta.st)/(exp(delta.st) + 
        1)^2 * ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) * 
        (log((-log(1 - (1 - p1)^teta))) * exp(delta.st)) + (-log(1 - 
        (1 - p2)^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - p2)^teta))) * exp(delta.st)))) - ((((-log(1 - (1 - 
    p1)^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - p1)^teta))) * exp(delta.st)) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p2)^teta))) * exp(delta.st)))) - 
    ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))) * 
        (exp(delta.st)/(exp(delta.st) + 1)^2))) * (log(((-log(1 - 
    (1 - p1)^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 
    1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)) + ((-log(1 - 
    (1 - p1)^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)) * (((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^teta))) * exp(delta.st)) + (-log(1 - 
    (1 - p2)^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    p2)^teta))) * exp(delta.st)))/((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1)) * (exp(delta.st)/(exp(delta.st) + 
    1)^2) + log(((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2 - exp(delta.st) * (2 * (exp(delta.st) * (exp(delta.st) + 
    1)))/((exp(delta.st) + 1)^2)^2)))) - exp(-((-log(1 - (1 - 
    p1)^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^teta))) * exp(delta.st)) + (-log(1 - 
    (1 - p2)^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    p2)^teta))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log(((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * (((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^teta))) * exp(delta.st)) + (-log(1 - 
    (1 - p2)^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    p2)^teta))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log(((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))))) - exp(-((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - 
    (1 - p1)^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^teta))) * exp(delta.st)) + (-log(1 - 
    (1 - p2)^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    p2)^teta))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log(((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))))) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^teta))) * exp(delta.st)) + (-log(1 - 
    (1 - p2)^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    p2)^teta))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (log(((-log(1 - (1 - p1)^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))))))




c.copula2.be1del <-(exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * ((((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - p1)^teta))^delta * 
    log((-log(1 - (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^delta * 
    log((-log(1 - (1 - p2)^teta))))) - ((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 1) * (log(((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)) * 
    (1/delta^2))) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
    1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - p1)^teta))))) + 
    ((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
        1) * log((-log(1 - (1 - p1)^teta))) * (delta * ((1 - 
        p1)^(teta - 1) * teta/(1 - (1 - p1)^teta))) + (-log(1 - 
        (1 - p1)^teta))^(delta - 1) * ((1 - p1)^(teta - 1) * 
        teta/(1 - (1 - p1)^teta))) - 1/delta^2 * ((-log(1 - (1 - 
        p1)^teta))^(delta - 1) * (delta * ((1 - p1)^(teta - 1) * 
        teta/(1 - (1 - p1)^teta)))))) - exp(-((-log(1 - (1 - 
    p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^delta * 
        log((-log(1 - (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^delta * 
        log((-log(1 - (1 - p2)^teta))))) - ((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta) * (log(((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)) * 
        (1/delta^2))) * (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
    (1 - p2)^teta))^delta)^((1/delta) - 1) * ((1/delta) * ((-log(1 - 
    (1 - p1)^teta))^(delta - 1) * (delta * ((1 - p1)^(teta - 
    1) * teta/(1 - (1 - p1)^teta)))))))) - exp(-((-log(1 - (1 - 
    p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^delta * log((-log(1 - 
    (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^delta * log((-log(1 - 
    (1 - p2)^teta))))) - ((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
    (1 - p2)^teta))^delta)^(1/delta) * (log(((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)) * (1/delta^2))))) * ((1/teta) * 
    (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^(delta - 
            1) * (delta * ((1 - p1)^(teta - 1) * teta/(1 - (1 - 
            p1)^teta)))))))))*exp(delta.st)


        




        
c.copula2.be2del <- (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * ((((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - p1)^teta))^delta * 
    log((-log(1 - (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^delta * 
    log((-log(1 - (1 - p2)^teta))))) - ((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 1) * (log(((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)) * 
    (1/delta^2))) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
    1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - p2)^teta))))) + 
    ((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
        1) * log((-log(1 - (1 - p2)^teta))) * (delta * ((1 - 
        p2)^(teta - 1) * teta/(1 - (1 - p2)^teta))) + (-log(1 - 
        (1 - p2)^teta))^(delta - 1) * ((1 - p2)^(teta - 1) * 
        teta/(1 - (1 - p2)^teta))) - 1/delta^2 * ((-log(1 - (1 - 
        p2)^teta))^(delta - 1) * (delta * ((1 - p2)^(teta - 1) * 
        teta/(1 - (1 - p2)^teta)))))) - exp(-((-log(1 - (1 - 
    p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^delta * 
        log((-log(1 - (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^delta * 
        log((-log(1 - (1 - p2)^teta))))) - ((-log(1 - (1 - p1)^teta))^delta + 
        (-log(1 - (1 - p2)^teta))^delta)^(1/delta) * (log(((-log(1 - 
        (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)) * 
        (1/delta^2))) * (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
    (1 - p2)^teta))^delta)^((1/delta) - 1) * ((1/delta) * ((-log(1 - 
    (1 - p2)^teta))^(delta - 1) * (delta * ((1 - p2)^(teta - 
    1) * teta/(1 - (1 - p2)^teta)))))))) - exp(-((-log(1 - (1 - 
    p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - p1)^teta))^delta * log((-log(1 - 
    (1 - p1)^teta))) + (-log(1 - (1 - p2)^teta))^delta * log((-log(1 - 
    (1 - p2)^teta))))) - ((-log(1 - (1 - p1)^teta))^delta + (-log(1 - 
    (1 - p2)^teta))^delta)^(1/delta) * (log(((-log(1 - (1 - p1)^teta))^delta + 
    (-log(1 - (1 - p2)^teta))^delta)) * (1/delta^2))))) * ((1/teta) * 
    (exp(-((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - p1)^teta))^delta + (-log(1 - (1 - p2)^teta))^delta)^((1/delta) - 
            1) * ((1/delta) * ((-log(1 - (1 - p2)^teta))^(delta - 
            1) * (delta * ((1 - p2)^(teta - 1) * teta/(1 - (1 - 
            p2)^teta)))))))))*exp(delta.st)




bit1.thdel <--(exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^((1/(exp(teta.st) + 1)) - 1) * 
    ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        ((((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^(((1/(exp(delta.st) + 1)) - 1) - 1) * (((1/(exp(delta.st) + 
            1)) - 1) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 1)))) * 
            exp(delta.st)) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
            1) * (log(((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * 
            ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
                1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
                1) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
                p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
                1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 
                1) - 1) * ((exp(delta.st) + 1) * ((1 - p2)^(exp(teta.st) + 
                1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - 
                p2)^(exp(teta.st) + 1)))))) + ((-log(1 - (1 - 
            p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
            (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
            1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - 
            (1 - p1)^(exp(teta.st) + 1)))^((exp(delta.st) + 1) - 
            1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 1)))) * 
            exp(delta.st)) * ((exp(delta.st) + 1) * ((1 - p1)^(exp(teta.st) + 
            1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - p1)^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * (exp(delta.st) * ((1 - p1)^(exp(teta.st) + 
            1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))) + ((-log(1 - (1 - p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))) * exp(delta.st)) * ((exp(delta.st) + 1) * ((1 - 
            p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
            (1 - p2)^(exp(teta.st) + 1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * (exp(delta.st) * 
            ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
                (1 - p2)^(exp(teta.st) + 1)))))) - exp(delta.st)/(exp(delta.st) + 
            1)^2 * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - p1)^(exp(teta.st) + 
            1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - p2)^(exp(teta.st) + 
            1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1))))))) - exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 1)))) * 
            exp(delta.st)) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
            (log(((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * 
        (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - p1)^(exp(teta.st) + 
            1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - p2)^(exp(teta.st) + 
            1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1))))))))) - exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 1)))) * exp(delta.st)) + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 1)))) * 
        exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))))) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - 
    p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^((exp(delta.st) + 
    1) - 1) * ((exp(delta.st) + 1) * ((1 - p1)^(exp(teta.st) + 
    1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
    1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 
    1) - 1) * ((exp(delta.st) + 1) * ((1 - p2)^(exp(teta.st) + 
    1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
    1))))))))) - (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^(1/(exp(teta.st) + 1)) * (exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 1)))) * exp(delta.st)) + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 1)))) * 
        exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)))/exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^((1/(exp(teta.st) + 1)) - 1) * 
    ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 1)))) * 
            exp(delta.st)) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
            (log(((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))))) * 
    (log(exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2))))


}


if(BivD=="BB6.90"){
   
  c.copula.be1 <- exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p1))^-teta))))))))


 
  c.copula.be2 <- 1 - exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
    (1 - p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - 
        p2)^-teta))))))))

  c.copula.theta <- (-(exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^(-1/teta) * (log(exp(-((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))) * 
    (1/teta^2)) - exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
    1) * (delta * ((1 - (1 - p1))^-teta * log((1 - (1 - p1)))/(1 - 
    (1 - (1 - p1))^-teta))) + (-log(1 - (1 - p2)^-teta))^-(delta + 
    1) * (delta * ((1 - p2)^-teta * log((1 - p2))/(1 - (1 - p2)^-teta))))))))))*(-exp(teta.st))


  c.copula.delta <- (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta) * 
        (log(((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)) * (1/delta^2)) - ((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-delta * 
        log((-log(1 - (1 - p2)^-teta))) + (-log(1 - (1 - (1 - 
        p1))^-teta))^-delta * log((-log(1 - (1 - (1 - p1))^-teta)))))))))*(-exp(delta.st))


 

 c.copula2.be1 <-  exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(((-1/delta) - 
        1) - 1) * (((-1/delta) - 1) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p1))^-teta))))) * ((-1/delta) * ((-log(1 - 
        (1 - (1 - p1))^-teta))^-(delta + 1) * (delta * ((1 - 
        (1 - p1))^-(teta + 1) * teta/(1 - (1 - (1 - p1))^-teta))))) + 
        ((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1 + 1) * 
                ((delta + 1) * ((1 - (1 - p1))^-(teta + 1) * 
                  teta/(1 - (1 - (1 - p1))^-teta))) * (delta * 
                ((1 - (1 - p1))^-(teta + 1) * teta/(1 - (1 - 
                  (1 - p1))^-teta))) - (-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
                1) * (delta * ((1 - (1 - p1))^-(teta + 1 + 1) * 
                (teta + 1) * teta/(1 - (1 - (1 - p1))^-teta) + 
                (1 - (1 - p1))^-(teta + 1) * teta * ((1 - (1 - 
                  p1))^-(teta + 1) * teta)/(1 - (1 - (1 - p1))^-teta)^2))))) - 
        exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
            (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
            1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
            (1 - (1 - p1))^-teta)))))) * (((-log(1 - (1 - (1 - 
            p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
            1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
            (1 - (1 - p1))^-teta)))))))) - exp(-((-log(1 - (1 - 
    (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
    1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - (1 - 
    (1 - p1))^-teta)))))))) * ((-1/teta) * (exp(-((-log(1 - (1 - 
    (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
        p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1) * (delta * 
            ((1 - (1 - p1))^-(teta + 1) * teta/(1 - (1 - (1 - 
                p1))^-teta))))))))
                  
 c.copula2.be2 <- -(exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
    (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - (1 - 
    p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
    1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta)))))))) * 
    ((-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - 
        p2)^-teta)))))))) + exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
    1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta)))))) * 
    (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
        p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - p2)^-teta))^-(delta + 1) * (delta * ((1 - 
            p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta)))))) + 
    exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * (delta * ((1 - p2)^-(teta + 1 + 1) * (teta + 1) * 
        teta/(1 - (1 - p2)^-teta) + (1 - p2)^-(teta + 1) * teta * 
        ((1 - p2)^-(teta + 1) * teta)/(1 - (1 - p2)^-teta)^2)) - 
        (-log(1 - (1 - p2)^-teta))^-(delta + 1 + 1) * ((delta + 
            1) * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta))) * 
            (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta))))) - 
        ((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^(((-1/delta) - 1) - 1) * 
            (((-1/delta) - 1) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
                1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - 
                (1 - p2)^-teta))))) * ((-1/delta) * ((-log(1 - 
            (1 - p2)^-teta))^-(delta + 1) * (delta * ((1 - p2)^-(teta + 
            1) * teta/(1 - (1 - p2)^-teta)))))))))
                  
                  
                 


c.copula2.be1be2 <- exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
    (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - (1 - 
    p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
    1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta)))))))) * 
    ((-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p1))^-teta)))))))) + exp(-((-log(1 - (1 - (1 - 
    p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
    1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta)))))) * 
    (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
        p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1) * (delta * 
            ((1 - (1 - p1))^-(teta + 1) * teta/(1 - (1 - (1 - 
                p1))^-teta)))))) - exp(-((-log(1 - (1 - (1 - 
    p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
        p2)^-teta))^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
        1) * ((-log(1 - (1 - p2)^-teta))^-(delta + 1) * (delta * 
        ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta))))) * 
        ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
            1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
            (1 - (1 - p1))^-teta))))))))




c.copula2.be1th <-((exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * (log(exp(-((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))) * 
    (1/teta^2)) - exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
    1) * (delta * ((1 - (1 - p1))^-teta * log((1 - (1 - p1)))/(1 - 
    (1 - (1 - p1))^-teta))) + (-log(1 - (1 - p2)^-teta))^-(delta + 
    1) * (delta * ((1 - p2)^-teta * log((1 - p2))/(1 - (1 - p2)^-teta))))))))) * 
    ((-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p1))^-teta)))))))) + exp(-((-log(1 - (1 - (1 - 
    p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 
    1) * (1/teta^2 * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
    1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - (1 - 
    (1 - p1))^-teta))))))) + (-1/teta) * (exp(-((-log(1 - (1 - 
    (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
        p2)^-teta))^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
        1) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1) * 
        (delta * ((1 - (1 - p1))^-teta * log((1 - (1 - p1)))/(1 - 
            (1 - (1 - p1))^-teta))) + (-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * (delta * ((1 - p2)^-teta * log((1 - p2))/(1 - (1 - 
        p2)^-teta))))) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p1))^-teta))))) + ((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * 
        ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
            1 + 1) * ((delta + 1) * ((1 - (1 - p1))^-teta * log((1 - 
            (1 - p1)))/(1 - (1 - (1 - p1))^-teta))) * (delta * 
            ((1 - (1 - p1))^-(teta + 1) * teta/(1 - (1 - (1 - 
                p1))^-teta))) + (-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
            1) * (delta * (((1 - (1 - p1))^-(teta + 1) - (1 - 
            (1 - p1))^-(teta + 1) * log((1 - (1 - p1))) * teta)/(1 - 
            (1 - (1 - p1))^-teta) - (1 - (1 - p1))^-(teta + 1) * 
            teta * ((1 - (1 - p1))^-teta * log((1 - (1 - p1))))/(1 - 
            (1 - (1 - p1))^-teta)^2))))) - exp(-((-log(1 - (1 - 
    (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
        p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1) * (delta * 
            ((1 - (1 - p1))^-teta * log((1 - (1 - p1)))/(1 - 
                (1 - (1 - p1))^-teta))) + (-log(1 - (1 - p2)^-teta))^-(delta + 
            1) * (delta * ((1 - p2)^-teta * log((1 - p2))/(1 - 
            (1 - p2)^-teta)))))) * (((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
    ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1) * (delta * 
        ((1 - (1 - p1))^-(teta + 1) * teta/(1 - (1 - (1 - p1))^-teta)))))))))*(-exp(teta.st))




c.copula2.be2th <-(-((exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
    (1 - p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * (log(exp(-((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))) * 
    (1/teta^2)) - exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
    1) * (delta * ((1 - (1 - p1))^-teta * log((1 - (1 - p1)))/(1 - 
    (1 - (1 - p1))^-teta))) + (-log(1 - (1 - p2)^-teta))^-(delta + 
    1) * (delta * ((1 - p2)^-teta * log((1 - p2))/(1 - (1 - p2)^-teta))))))))) * 
    ((-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - 
        p2)^-teta)))))))) + exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 
    1) * (1/teta^2 * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
    1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta))))))) + 
    (-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(((-1/delta) - 
        1) - 1) * (((-1/delta) - 1) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^-teta * log((1 - (1 - p1)))/(1 - 
        (1 - (1 - p1))^-teta))) + (-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * (delta * ((1 - p2)^-teta * log((1 - p2))/(1 - (1 - 
        p2)^-teta))))) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - 
        p2)^-teta))))) + ((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * 
        ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 1 + 
            1) * ((delta + 1) * ((1 - p2)^-teta * log((1 - p2))/(1 - 
            (1 - p2)^-teta))) * (delta * ((1 - p2)^-(teta + 1) * 
            teta/(1 - (1 - p2)^-teta))) + (-log(1 - (1 - p2)^-teta))^-(delta + 
            1) * (delta * (((1 - p2)^-(teta + 1) - (1 - p2)^-(teta + 
            1) * log((1 - p2)) * teta)/(1 - (1 - p2)^-teta) - 
            (1 - p2)^-(teta + 1) * teta * ((1 - p2)^-teta * log((1 - 
                p2)))/(1 - (1 - p2)^-teta)^2))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1) * 
                (delta * ((1 - (1 - p1))^-teta * log((1 - (1 - 
                  p1)))/(1 - (1 - (1 - p1))^-teta))) + (-log(1 - 
                (1 - p2)^-teta))^-(delta + 1) * (delta * ((1 - 
                p2)^-teta * log((1 - p2))/(1 - (1 - p2)^-teta)))))) * 
        (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - p2)^-teta))^-(delta + 1) * (delta * 
                ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta))))))))))*(-exp(teta.st))







bit1.th2 <--((exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - 
    p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - 
    p2)^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
    ((-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-(delta + 1) * (delta * 
        ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
            (1 - p2)^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))^-(delta + 1) * (delta * ((1 - 
        (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * 
        exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 1))))))))) - 
    exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + 
        (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))^((1/(exp(teta.st) + 
        1)) - 1) * (log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))) * 
        (exp(teta.st)/(exp(teta.st) + 1)^2))) * ((1/(exp(teta.st) + 
    1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + 
        (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
        (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1))))))))) + exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^-delta)^(-1/delta)) * (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-(delta + 
    1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
    exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 1)))) + (-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^-(delta + 1) * (delta * 
    ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * 
        exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 1))))))) * 
    (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + 
        (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
        (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1))))))) + exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + 
        (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * (((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st)) * (log((1 - p2)) * 
        exp(teta.st)) + (1 - p2)^(exp(teta.st) + 1) * (log((1 - 
        p2)) * exp(teta.st)))/(1 - (1 - p2)^(exp(teta.st) + 1)) + 
        (1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)) * 
            ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)))/(1 - 
            (1 - p2)^(exp(teta.st) + 1))^2)) - (-log(1 - (1 - 
        p2)^(exp(teta.st) + 1)))^-(delta + 1 + 1) * ((delta + 
        1) * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
        exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 1)))) * (delta * 
        ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
            (1 - p2)^(exp(teta.st) + 1)))) + ((-log(1 - (1 - 
        (1 - p1))^(exp(teta.st) + 1)))^-(delta + 1) * (delta * 
        (((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - 
            p1))) * exp(teta.st)) * (log((1 - (1 - p1))) * exp(teta.st)) + 
            (1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - 
                p1))) * exp(teta.st)))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)) + (1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
            (1 - p1))) * exp(teta.st)) * ((1 - (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (1 - p1))) * exp(teta.st)))/(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1))^2)) - (-log(1 - (1 - 
        (1 - p1))^(exp(teta.st) + 1)))^-(delta + 1 + 1) * ((delta + 
        1) * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - 
        p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 1) * 
        (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1))))))) - ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(((-1/delta) - 
        1) - 1) * (((-1/delta) - 1) * ((-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
        (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))))) * ((-1/delta) * ((-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
        (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))))))) - exp(teta.st)/(exp(teta.st) + 1)^2 * (exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - 
    p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - 
    p2)^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
    ((-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-(delta + 1) * (delta * 
        ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
            (1 - p2)^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))^-(delta + 1) * (delta * ((1 - 
        (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * 
        exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 1))))))))) - 
    ((exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + 
        (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))^((1/(exp(teta.st) + 
        1)) - 1) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - 
        (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 1) * 
        ((-1/delta) * ((-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-(delta + 
            1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
            p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^-(delta + 1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1))))))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))^(1/(exp(teta.st) + 
        1)) * (log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))) * 
        (exp(teta.st)/(exp(teta.st) + 1)^2))) * (log(exp(-((-log(1 - 
        (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))) * 
        (exp(teta.st)/(exp(teta.st) + 1)^2)) + exp(-((-log(1 - 
        (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))^(1/(exp(teta.st) + 
        1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + 
            (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^-(delta + 1) * (delta * ((1 - p2)^(exp(teta.st) + 
            1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^-(delta + 1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))))))/exp(-((-log(1 - 
        (1 - (1 - p1))^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * (exp(teta.st)/(exp(teta.st) + 
        1)^2) + log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^-delta + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^-delta)^(-1/delta))) * 
        (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * (2 * 
            (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
            1)^2)^2))))




bit1.del2 <-exp(-((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^((-1/teta) - 1) * ((-1/teta) * (exp(-((-log(1 - (1 - 
    (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * ((((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 
    1)) - 1) - 1) * (((1/(exp(delta.st) + 1)) - 1) * ((-log(1 - 
    (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - (1 - p1))^-teta))) * exp(delta.st)) + (-log(1 - (1 - 
    p2)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p2)^-teta))) * 
    exp(delta.st)))) - ((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * (log(((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - 
    p1))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
    p1))^-teta))) * exp(delta.st)) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p2)^-teta))) * exp(delta.st)))) + 
    ((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - p2)^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - 
        (1 - p1))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p1))^-teta))) * exp(delta.st)) * (log((-log(1 - 
        (1 - (1 - p1))^-teta))) * exp(delta.st)) + (-log(1 - 
        (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p1))^-teta))) * exp(delta.st)) + ((-log(1 - 
        (1 - p2)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - p2)^-teta))) * exp(delta.st)) * (log((-log(1 - (1 - 
        p2)^-teta))) * exp(delta.st)) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 
        1) * (log((-log(1 - (1 - p2)^-teta))) * exp(delta.st)))) - 
        exp(delta.st)/(exp(delta.st) + 1)^2 * ((-log(1 - (1 - 
            (1 - p1))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
            (1 - (1 - p1))^-teta))) * exp(delta.st)) + (-log(1 - 
            (1 - p2)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
            (1 - p2)^-teta))) * exp(delta.st)))) - ((((-log(1 - 
    (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
    p2)^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
    1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p1))^-teta))) * exp(delta.st)) + 
    (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - p2)^-teta))) * exp(delta.st)))) - ((-log(1 - (1 - 
    (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * (log(((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2)) + ((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)) * (((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p1))^-teta))) * exp(delta.st)) + 
    (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - p2)^-teta))) * exp(delta.st)))/((-log(1 - (1 - (1 - 
    p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 
    1)) * (exp(delta.st)/(exp(delta.st) + 1)^2) + log(((-log(1 - 
    (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
    p2)^-teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2 - exp(delta.st) * (2 * (exp(delta.st) * (exp(delta.st) + 
    1)))/((exp(delta.st) + 1)^2)^2)))) - exp(-((-log(1 - (1 - 
    (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - 
    p1))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
    p1))^-teta))) * exp(delta.st)) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p2)^-teta))) * exp(delta.st)))) - 
    ((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - p2)^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))) * 
        (exp(delta.st)/(exp(delta.st) + 1)^2))) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
    p2)^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
    1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p1))^-teta))) * exp(delta.st)) + 
    (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - p2)^-teta))) * exp(delta.st)))) - ((-log(1 - (1 - 
    (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))))) - exp(-((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^(((-1/teta) - 1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - 
    (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
    p2)^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
    (((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 1) + 
        (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - 
        (1 - p1))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p1))^-teta))) * exp(delta.st)) + (-log(1 - 
        (1 - p2)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - p2)^-teta))) * exp(delta.st)))) - ((-log(1 - (1 - 
        (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
        p2)^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))) * 
        (exp(delta.st)/(exp(delta.st) + 1)^2))))) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1))) * (((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - 
        (1 - p1))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p1))^-teta))) * exp(delta.st)) + (-log(1 - 
        (1 - p2)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - p2)^-teta))) * exp(delta.st)))) - ((-log(1 - (1 - 
        (1 - p1))^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
        p2)^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - (1 - p1))^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^-teta))^(exp(delta.st) + 1))) * 
        (exp(delta.st)/(exp(delta.st) + 1)^2)))))




c.copula2.be1del <-(exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta)) * ((((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * (log(((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)) * (1/delta^2)) - 
        ((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^(((-1/delta) - 1) - 1) * 
            (((-1/delta) - 1) * ((-log(1 - (1 - p2)^-teta))^-delta * 
                log((-log(1 - (1 - p2)^-teta))) + (-log(1 - (1 - 
                (1 - p1))^-teta))^-delta * log((-log(1 - (1 - 
                (1 - p1))^-teta)))))) * ((-1/delta) * ((-log(1 - 
        (1 - (1 - p1))^-teta))^-(delta + 1) * (delta * ((1 - 
        (1 - p1))^-(teta + 1) * teta/(1 - (1 - (1 - p1))^-teta))))) + 
        ((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * (1/delta^2 * 
            ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1) * 
                (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
                  (1 - (1 - p1))^-teta)))) + (-1/delta) * ((-log(1 - 
            (1 - (1 - p1))^-teta))^-(delta + 1) * ((1 - (1 - 
            p1))^-(teta + 1) * teta/(1 - (1 - (1 - p1))^-teta)) - 
            (-log(1 - (1 - (1 - p1))^-teta))^-(delta + 1) * log((-log(1 - 
                (1 - (1 - p1))^-teta))) * (delta * ((1 - (1 - 
                p1))^-(teta + 1) * teta/(1 - (1 - (1 - p1))^-teta)))))) - 
        exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
            (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta) * 
            (log(((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
                (1 - p2)^-teta))^-delta)) * (1/delta^2)) - ((-log(1 - 
            (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-delta * 
            log((-log(1 - (1 - p2)^-teta))) + (-log(1 - (1 - 
            (1 - p1))^-teta))^-delta * log((-log(1 - (1 - (1 - 
            p1))^-teta)))))) * (((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
            (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
            1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
            (1 - (1 - p1))^-teta)))))))) - exp(-((-log(1 - (1 - 
    (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
    (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta) * 
    (log(((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)) * (1/delta^2)) - ((-log(1 - 
    (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-delta * log((-log(1 - 
    (1 - p2)^-teta))) + (-log(1 - (1 - (1 - p1))^-teta))^-delta * 
    log((-log(1 - (1 - (1 - p1))^-teta)))))))) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p1))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p1))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p1))^-teta))))))))
)*(-exp(delta.st))


        



        
c.copula2.be2del <- (-(exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - 
    p2)^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta)) * ((((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * (log(((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)) * (1/delta^2)) - 
        ((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^(((-1/delta) - 1) - 1) * 
            (((-1/delta) - 1) * ((-log(1 - (1 - p2)^-teta))^-delta * 
                log((-log(1 - (1 - p2)^-teta))) + (-log(1 - (1 - 
                (1 - p1))^-teta))^-delta * log((-log(1 - (1 - 
                (1 - p1))^-teta)))))) * ((-1/delta) * ((-log(1 - 
        (1 - p2)^-teta))^-(delta + 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta/(1 - (1 - p2)^-teta))))) + ((-log(1 - (1 - 
        (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * (1/delta^2 * ((-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - 
        p2)^-teta)))) + (-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 
        1) * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta)) - 
        (-log(1 - (1 - p2)^-teta))^-(delta + 1) * log((-log(1 - 
            (1 - p2)^-teta))) * (delta * ((1 - p2)^-(teta + 1) * 
            teta/(1 - (1 - p2)^-teta)))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)^(-1/delta) * (log(((-log(1 - 
            (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)) * 
            (1/delta^2)) - ((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
            (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-delta * 
            log((-log(1 - (1 - p2)^-teta))) + (-log(1 - (1 - 
            (1 - p1))^-teta))^-delta * log((-log(1 - (1 - (1 - 
            p1))^-teta)))))) * (((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * 
        ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-(delta + 1) * 
            (delta * ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta)))))))) - 
    exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
        (1 - p2)^-teta))^-delta)^(-1/delta))^(((-1/teta) - 1) - 
        1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
        (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta) * 
        (log(((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
            (1 - p2)^-teta))^-delta)) * (1/delta^2)) - ((-log(1 - 
        (1 - (1 - p1))^-teta))^-delta + (-log(1 - (1 - p2)^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p2)^-teta))^-delta * 
        log((-log(1 - (1 - p2)^-teta))) + (-log(1 - (1 - (1 - 
        p1))^-teta))^-delta * log((-log(1 - (1 - (1 - p1))^-teta)))))))) * 
        ((-1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^-teta))^-delta + 
            (-log(1 - (1 - p2)^-teta))^-delta)^(-1/delta)) * 
            (((-log(1 - (1 - (1 - p1))^-teta))^-delta + (-log(1 - 
                (1 - p2)^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
                ((-log(1 - (1 - p2)^-teta))^-(delta + 1) * (delta * 
                  ((1 - p2)^-(teta + 1) * teta/(1 - (1 - p2)^-teta)))))))))
)*(-exp(delta.st))




bit1.thdel <-exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^((1/(exp(teta.st) + 1)) - 1) * 
    ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        ((((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^(((1/(exp(delta.st) + 1)) - 1) - 1) * (((1/(exp(delta.st) + 
            1)) - 1) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
            p1))^(exp(teta.st) + 1)))) * exp(delta.st)) + (-log(1 - 
            (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) * 
            (log((-log(1 - (1 - p2)^(exp(teta.st) + 1)))) * exp(delta.st)))) - 
            ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))^((1/(exp(delta.st) + 1)) - 1) * (log(((-log(1 - 
                (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * 
            ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
                1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
                1) * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
                (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
                1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 
                1) - 1) * ((exp(delta.st) + 1) * ((1 - p2)^(exp(teta.st) + 
                1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - 
                p2)^(exp(teta.st) + 1)))))) + ((-log(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
            (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * (log((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) * exp(delta.st)) * ((exp(delta.st) + 1) * ((1 - 
            (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * 
            exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * (exp(delta.st) * 
            ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - 
                p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
                1)))) + ((-log(1 - (1 - p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))) * exp(delta.st)) * ((exp(delta.st) + 1) * ((1 - 
            p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
            (1 - p2)^(exp(teta.st) + 1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * (exp(delta.st) * 
            ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
                (1 - p2)^(exp(teta.st) + 1)))))) - exp(delta.st)/(exp(delta.st) + 
            1)^2 * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
            1) * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
            (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - p2)^(exp(teta.st) + 
            1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - p2)^(exp(teta.st) + 
            1))))))) - exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) * exp(delta.st)) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))) * exp(delta.st)))) - ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
            (log(((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * 
        (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - 
            p2)^(exp(teta.st) + 1)))^((exp(delta.st) + 1) - 1) * 
            ((exp(delta.st) + 1) * ((1 - p2)^(exp(teta.st) + 
                1) * (log((1 - p2)) * exp(teta.st))/(1 - (1 - 
                p2)^(exp(teta.st) + 1))))))))) - exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)))^(((1/(exp(teta.st) + 1)) - 
    1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))) * exp(delta.st)) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))) * exp(delta.st)))) - ((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))))) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^((exp(delta.st) + 
    1) - 1) * ((exp(delta.st) + 1) * ((1 - (1 - p1))^(exp(teta.st) + 
    1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - p2)^(exp(teta.st) + 
    1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * ((1 - 
    p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))/(1 - 
    (1 - p2)^(exp(teta.st) + 1))))))))) - (exp(-((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^(1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
    1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))) * 
    exp(delta.st)) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 1)))) * exp(delta.st)))) - 
    ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 1))) * 
        (exp(delta.st)/(exp(delta.st) + 1)^2)))/exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^((1/(exp(teta.st) + 1)) - 1) * 
    ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) * exp(delta.st)) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))) * exp(delta.st)))) - ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) + (-log(1 - (1 - p2)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
            (log(((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))))) * 
    (log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1) + (-log(1 - (1 - p2)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)))


}


if(BivD=="BB6.180"){

   
  c.copula.be1 <- 1 - exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p1))^(teta - 1) * teta/(1 - 
        (1 - (1 - p1))^teta))))))))

 
  c.copula.be2 <- 1 - exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p2))^(teta - 1) * teta/(1 - 
        (1 - (1 - p2))^teta))))))))



  c.copula.theta <-(-(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^(1/teta) * (log(exp(-((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))) * 
    (1/teta^2)) + exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
    (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 
    1) * ((1/teta) * (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
    (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
    1) * (delta * ((1 - (1 - p1))^teta * log((1 - (1 - p1)))/(1 - 
    (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - p2))^teta))^(delta - 
    1) * (delta * ((1 - (1 - p2))^teta * log((1 - (1 - p2)))/(1 - 
    (1 - (1 - p2))^teta))))))))))*exp(teta.st)


  c.copula.delta <- (-(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^delta * 
        log((-log(1 - (1 - (1 - p1))^teta))) + (-log(1 - (1 - 
        (1 - p2))^teta))^delta * log((-log(1 - (1 - (1 - p2))^teta))))) - 
        ((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^(1/delta) * (log(((-log(1 - 
            (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - 
            p2))^teta))^delta)) * (1/delta^2)))))))*exp(delta.st)




 c.copula2.be1 <- -(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p1))^(teta - 1) * teta/(1 - 
        (1 - (1 - p1))^teta))))) * ((1/delta) * ((-log(1 - (1 - 
        (1 - p1))^teta))^(delta - 1) * (delta * ((1 - (1 - p1))^(teta - 
        1) * teta/(1 - (1 - (1 - p1))^teta))))) + ((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^((delta - 
        1) - 1) * ((delta - 1) * ((1 - (1 - p1))^(teta - 1) * 
        teta/(1 - (1 - (1 - p1))^teta))) * (delta * ((1 - (1 - 
        p1))^(teta - 1) * teta/(1 - (1 - (1 - p1))^teta))) + 
        (-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p1))^((teta - 1) - 1) * (teta - 1) * teta/(1 - 
                (1 - (1 - p1))^teta) + (1 - (1 - p1))^(teta - 
                1) * teta * ((1 - (1 - p1))^(teta - 1) * teta)/(1 - 
                (1 - (1 - p1))^teta)^2))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
            ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
                ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - 
                  p1))^teta)))))) * (((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
            1) * (delta * ((1 - (1 - p1))^(teta - 1) * teta/(1 - 
            (1 - (1 - p1))^teta)))))))) - exp(-((-log(1 - (1 - 
    (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
    (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
    1) * (delta * ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - 
    (1 - p1))^teta)))))))) * ((1/teta) * (exp(-((-log(1 - (1 - 
    (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - 
                p1))^teta)))))))))


                  
 c.copula2.be2 <- -(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p2))^(teta - 1) * teta/(1 - 
        (1 - (1 - p2))^teta))))) * ((1/delta) * ((-log(1 - (1 - 
        (1 - p2))^teta))^(delta - 1) * (delta * ((1 - (1 - p2))^(teta - 
        1) * teta/(1 - (1 - (1 - p2))^teta))))) + ((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^((delta - 
        1) - 1) * ((delta - 1) * ((1 - (1 - p2))^(teta - 1) * 
        teta/(1 - (1 - (1 - p2))^teta))) * (delta * ((1 - (1 - 
        p2))^(teta - 1) * teta/(1 - (1 - (1 - p2))^teta))) + 
        (-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p2))^((teta - 1) - 1) * (teta - 1) * teta/(1 - 
                (1 - (1 - p2))^teta) + (1 - (1 - p2))^(teta - 
                1) * teta * ((1 - (1 - p2))^(teta - 1) * teta)/(1 - 
                (1 - (1 - p2))^teta)^2))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
            ((-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
                ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - 
                  p2))^teta)))))) * (((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
            1) * (delta * ((1 - (1 - p2))^(teta - 1) * teta/(1 - 
            (1 - (1 - p2))^teta)))))))) - exp(-((-log(1 - (1 - 
    (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
    (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
    1) * (delta * ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - 
    (1 - p2))^teta)))))))) * ((1/teta) * (exp(-((-log(1 - (1 - 
    (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - 
                p2))^teta)))))))))

                  
                  




c.copula2.be1be2 <- -(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p2))^(teta - 1) * teta/(1 - 
        (1 - (1 - p2))^teta))))) * ((1/delta) * ((-log(1 - (1 - 
        (1 - p1))^teta))^(delta - 1) * (delta * ((1 - (1 - p1))^(teta - 
        1) * teta/(1 - (1 - (1 - p1))^teta)))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
            ((-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
                ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - 
                  p2))^teta)))))) * (((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
            1) * (delta * ((1 - (1 - p1))^(teta - 1) * teta/(1 - 
            (1 - (1 - p1))^teta)))))))) - exp(-((-log(1 - (1 - 
    (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
    (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
    1) * (delta * ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - 
    (1 - p2))^teta)))))))) * ((1/teta) * (exp(-((-log(1 - (1 - 
    (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - 
                p1))^teta)))))))))





c.copula2.be1th <-(-(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p1))^teta * log((1 - (1 - p1)))/(1 - 
        (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - p2))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p2))^teta * log((1 - (1 - p2)))/(1 - 
        (1 - (1 - p2))^teta))))) * ((1/delta) * ((-log(1 - (1 - 
        (1 - p1))^teta))^(delta - 1) * (delta * ((1 - (1 - p1))^(teta - 
        1) * teta/(1 - (1 - (1 - p1))^teta))))) + ((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^((delta - 
        1) - 1) * ((delta - 1) * ((1 - (1 - p1))^teta * log((1 - 
        (1 - p1)))/(1 - (1 - (1 - p1))^teta))) * (delta * ((1 - 
        (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - p1))^teta))) + 
        (-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
            (((1 - (1 - p1))^(teta - 1) * log((1 - (1 - p1))) * 
                teta + (1 - (1 - p1))^(teta - 1))/(1 - (1 - (1 - 
                p1))^teta) + (1 - (1 - p1))^(teta - 1) * teta * 
                ((1 - (1 - p1))^teta * log((1 - (1 - p1))))/(1 - 
                (1 - (1 - p1))^teta)^2))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
            ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
                ((1 - (1 - p1))^teta * log((1 - (1 - p1)))/(1 - 
                  (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - 
                p2))^teta))^(delta - 1) * (delta * ((1 - (1 - 
                p2))^teta * log((1 - (1 - p2)))/(1 - (1 - (1 - 
                p2))^teta)))))) * (((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
            1) * (delta * ((1 - (1 - p1))^(teta - 1) * teta/(1 - 
            (1 - (1 - p1))^teta))))))) - 1/teta^2 * (exp(-((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - 
                p1))^teta)))))))) - (exp(-((-log(1 - (1 - (1 - 
    p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 
    1) * (log(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
    (1 - (1 - p2))^teta))^delta)^(1/delta))) * (1/teta^2)) + 
    exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^(1/delta))^(((1/teta) - 1) - 
        1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p1))^teta * log((1 - (1 - p1)))/(1 - 
        (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - p2))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p2))^teta * log((1 - (1 - p2)))/(1 - 
        (1 - (1 - p2))^teta))))))))) * ((1/teta) * (exp(-((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - 
                p1))^teta))))))))))*exp(teta.st)


        


c.copula2.be2th <-(-(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p1))^teta * log((1 - (1 - p1)))/(1 - 
        (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - p2))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p2))^teta * log((1 - (1 - p2)))/(1 - 
        (1 - (1 - p2))^teta))))) * ((1/delta) * ((-log(1 - (1 - 
        (1 - p2))^teta))^(delta - 1) * (delta * ((1 - (1 - p2))^(teta - 
        1) * teta/(1 - (1 - (1 - p2))^teta))))) + ((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^((delta - 
        1) - 1) * ((delta - 1) * ((1 - (1 - p2))^teta * log((1 - 
        (1 - p2)))/(1 - (1 - (1 - p2))^teta))) * (delta * ((1 - 
        (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - p2))^teta))) + 
        (-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
            (((1 - (1 - p2))^(teta - 1) * log((1 - (1 - p2))) * 
                teta + (1 - (1 - p2))^(teta - 1))/(1 - (1 - (1 - 
                p2))^teta) + (1 - (1 - p2))^(teta - 1) * teta * 
                ((1 - (1 - p2))^teta * log((1 - (1 - p2))))/(1 - 
                (1 - (1 - p2))^teta)^2))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
            ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
                ((1 - (1 - p1))^teta * log((1 - (1 - p1)))/(1 - 
                  (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - 
                p2))^teta))^(delta - 1) * (delta * ((1 - (1 - 
                p2))^teta * log((1 - (1 - p2)))/(1 - (1 - (1 - 
                p2))^teta)))))) * (((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
            1) * (delta * ((1 - (1 - p2))^(teta - 1) * teta/(1 - 
            (1 - (1 - p2))^teta))))))) - 1/teta^2 * (exp(-((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - 
                p2))^teta)))))))) - (exp(-((-log(1 - (1 - (1 - 
    p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 
    1) * (log(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
    (1 - (1 - p2))^teta))^delta)^(1/delta))) * (1/teta^2)) + 
    exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^(1/delta))^(((1/teta) - 1) - 
        1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p1))^teta * log((1 - (1 - p1)))/(1 - 
        (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - p2))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p2))^teta * log((1 - (1 - p2)))/(1 - 
        (1 - (1 - p2))^teta))))))))) * ((1/teta) * (exp(-((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - 
                p2))^teta))))))))))*exp(teta.st)







bit1.th2 <--(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta))^(1/(exp(teta.st) + 
    1)) * (log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * (2 * 
        (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
        1)^2)^2) - exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1)))))))/exp(-((-log(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^delta)^(1/delta)) * (exp(teta.st)/(exp(teta.st) + 1)^2)) - 
    (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + 
        (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta))^(1/(exp(teta.st) + 
        1)) * (log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^delta)^(1/delta))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)) + exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^delta)^(1/delta))^((1/(exp(teta.st) + 1)) - 1) * 
        ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - (1 - 
            p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
            (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + 
                (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^((1/delta) - 
                1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
                1)))^(delta - 1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 
                1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - 
                (1 - (1 - p1))^(exp(teta.st) + 1)))) + (-log(1 - 
                (1 - (1 - p2))^(exp(teta.st) + 1)))^(delta - 
                1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
                1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - 
                (1 - (1 - p2))^(exp(teta.st) + 1)))))))))) * 
        (log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + 
            (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta))) * 
            (exp(teta.st)/(exp(teta.st) + 1)^2)) + (exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
    (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
    (1 - p2))^(exp(teta.st) + 1)))^delta)^(((1/delta) - 1) - 
    1) * (((1/delta) - 1) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^(delta - 1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 
    1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(delta - 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
    1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))))) * ((1/delta) * ((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^(delta - 1) * (delta * ((1 - 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^(delta - 1) * (delta * ((1 - (1 - 
    p2))^(exp(teta.st) + 1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))))) + ((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^delta)^((1/delta) - 1) * ((1/delta) * 
    ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^((delta - 
        1) - 1) * ((delta - 1) * ((1 - (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * (((1 - (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (1 - p1))) * exp(teta.st)) * (log((1 - 
        (1 - p1))) * exp(teta.st)) + (1 - (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (1 - p1))) * exp(teta.st)))/(1 - (1 - 
        (1 - p1))^(exp(teta.st) + 1)) + (1 - (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (1 - p1))) * exp(teta.st)) * ((1 - (1 - 
        p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * exp(teta.st)))/(1 - 
        (1 - (1 - p1))^(exp(teta.st) + 1))^2)) + ((-log(1 - (1 - 
        (1 - p2))^(exp(teta.st) + 1)))^((delta - 1) - 1) * ((delta - 
        1) * ((1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - (1 - 
        p2))) * exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 1) * 
        (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(delta - 
        1) * (delta * (((1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - 
        (1 - p2))) * exp(teta.st)) * (log((1 - (1 - p2))) * exp(teta.st)) + 
        (1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - (1 - p2))) * 
            exp(teta.st)))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)) + (1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - (1 - 
        p2))) * exp(teta.st)) * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st)))/(1 - (1 - 
        (1 - p2))^(exp(teta.st) + 1))^2)))))) - exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
    (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
    (1 - p2))^(exp(teta.st) + 1)))^delta)^((1/delta) - 1) * ((1/delta) * 
    ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(delta - 1) * 
        (delta * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
            (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1))))))) * (((-log(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^delta)^((1/delta) - 1) * ((1/delta) * ((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^(delta - 1) * (delta * ((1 - 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^(delta - 1) * (delta * ((1 - (1 - 
    p2))^(exp(teta.st) + 1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))))))) - exp(teta.st)/(exp(teta.st) + 
    1)^2 * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1))))))))) - (exp(-((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^delta)^(1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * (log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) + exp(-((-log(1 - (1 - 
    (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^delta)^(1/delta))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
    (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - (1 - 
    (1 - p2))^(exp(teta.st) + 1)))^delta)^((1/delta) - 1) * ((1/delta) * 
    ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(delta - 1) * 
        (delta * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
            (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1)))))))))) * ((1/(exp(teta.st) + 
    1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^delta + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(delta - 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1)))))))))))




bit1.del2 <--(exp(-((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^((1/teta) - 1) * ((1/teta) * (exp(-((-log(1 - (1 - (1 - 
    p1))^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * ((((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 
    1)) - 1) - 1) * (((1/(exp(delta.st) + 1)) - 1) * ((-log(1 - 
    (1 - (1 - p1))^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - (1 - p1))^teta))) * exp(delta.st)) + (-log(1 - (1 - 
    (1 - p2))^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    (1 - p2))^teta))) * exp(delta.st)))) - ((-log(1 - (1 - (1 - 
    p1))^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
    1))^((1/(exp(delta.st) + 1)) - 1) * (log(((-log(1 - (1 - 
    (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - 
    p2))^teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
    1)^2))) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - 
    p1))^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
    p1))^teta))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p2))^teta))) * exp(delta.st)))) + 
    ((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - 
        (1 - p1))^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p1))^teta))) * exp(delta.st)) * (log((-log(1 - 
        (1 - (1 - p1))^teta))) * exp(delta.st)) + (-log(1 - (1 - 
        (1 - p1))^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p1))^teta))) * exp(delta.st)) + ((-log(1 - 
        (1 - (1 - p2))^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^teta))) * exp(delta.st)) * (log((-log(1 - 
        (1 - (1 - p2))^teta))) * exp(delta.st)) + (-log(1 - (1 - 
        (1 - p2))^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^teta))) * exp(delta.st)))) - exp(delta.st)/(exp(delta.st) + 
        1)^2 * ((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) * (log((-log(1 - (1 - (1 - p1))^teta))) * exp(delta.st)) + 
        (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1) * 
            (log((-log(1 - (1 - (1 - p2))^teta))) * exp(delta.st)))) - 
    ((((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + 
        (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - 
        (1 - p1))^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p1))^teta))) * exp(delta.st)) + (-log(1 - (1 - 
        (1 - p2))^teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^teta))) * exp(delta.st)))) - ((-log(1 - 
        (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * (log(((-log(1 - 
        (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2)) + ((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) * (log((-log(1 - (1 - (1 - p1))^teta))) * exp(delta.st)) + 
        (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1) * 
            (log((-log(1 - (1 - (1 - p2))^teta))) * exp(delta.st)))/((-log(1 - 
        (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^teta))^(exp(delta.st) + 1)) * (exp(delta.st)/(exp(delta.st) + 
        1)^2) + log(((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2 - exp(delta.st) * 
        (2 * (exp(delta.st) * (exp(delta.st) + 1)))/((exp(delta.st) + 
        1)^2)^2)))) - exp(-((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - 
    p1))^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
    p1))^teta))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p2))^teta))) * exp(delta.st)))) - 
    ((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * (((-log(1 - 
    (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
    (1 - p2))^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - 
    p1))^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
    p1))^teta))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p2))^teta))) * exp(delta.st)))) - 
    ((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))))) - exp(-((-log(1 - 
    (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
    (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - 
    (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
    (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - 
    p1))^teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
    p1))^teta))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p2))^teta))) * exp(delta.st)))) - 
    ((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))))) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + 
        (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1))) * (((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * ((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 
        1) * (log((-log(1 - (1 - (1 - p1))^teta))) * exp(delta.st)) + 
        (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 1) * 
            (log((-log(1 - (1 - (1 - p2))^teta))) * exp(delta.st)))) - 
        ((-log(1 - (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + 
            (-log(1 - (1 - (1 - p2))^teta))^(exp(delta.st) + 
                1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - 
            (1 - (1 - p1))^teta))^(exp(delta.st) + 1) + (-log(1 - 
            (1 - (1 - p2))^teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
            1)^2))))))


c.copula2.be1del <- (-(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * ((((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - (1 - p1))^teta))^delta * 
        log((-log(1 - (1 - (1 - p1))^teta))) + (-log(1 - (1 - 
        (1 - p2))^teta))^delta * log((-log(1 - (1 - (1 - p2))^teta))))) - 
        ((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * (log(((-log(1 - 
            (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - 
            p2))^teta))^delta)) * (1/delta^2))) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - 
                p1))^teta))))) + ((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
            1) * log((-log(1 - (1 - (1 - p1))^teta))) * (delta * 
            ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - 
                p1))^teta))) + (-log(1 - (1 - (1 - p1))^teta))^(delta - 
            1) * ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - 
            (1 - p1))^teta))) - 1/delta^2 * ((-log(1 - (1 - (1 - 
            p1))^teta))^(delta - 1) * (delta * ((1 - (1 - p1))^(teta - 
            1) * teta/(1 - (1 - (1 - p1))^teta)))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
            ((-log(1 - (1 - (1 - p1))^teta))^delta * log((-log(1 - 
                (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - 
                p2))^teta))^delta * log((-log(1 - (1 - (1 - p2))^teta))))) - 
            ((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
                (1 - (1 - p2))^teta))^delta)^(1/delta) * (log(((-log(1 - 
                (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
                (1 - p2))^teta))^delta)) * (1/delta^2))) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p1))^(teta - 1) * teta/(1 - 
        (1 - (1 - p1))^teta)))))))) - exp(-((-log(1 - (1 - (1 - 
    p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
    (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^delta * 
    log((-log(1 - (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - 
    p2))^teta))^delta * log((-log(1 - (1 - (1 - p2))^teta))))) - 
    ((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - 
        p2))^teta))^delta)^(1/delta) * (log(((-log(1 - (1 - (1 - 
        p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)) * 
        (1/delta^2))))) * ((1/teta) * (exp(-((-log(1 - (1 - (1 - 
    p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p1))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p1))^(teta - 1) * teta/(1 - (1 - (1 - 
                p1))^teta))))))))))*exp(delta.st)


        



        
c.copula2.be2del <- (-(exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
    (1 - p2))^teta))^delta)^(1/delta))^((1/teta) - 1) * ((1/teta) * 
    (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
        (1 - (1 - p2))^teta))^delta)^(1/delta)) * ((((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((-log(1 - (1 - (1 - p1))^teta))^delta * 
        log((-log(1 - (1 - (1 - p1))^teta))) + (-log(1 - (1 - 
        (1 - p2))^teta))^delta * log((-log(1 - (1 - (1 - p2))^teta))))) - 
        ((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * (log(((-log(1 - 
            (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - 
            p2))^teta))^delta)) * (1/delta^2))) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - 
                p2))^teta))))) + ((-log(1 - (1 - (1 - p1))^teta))^delta + 
        (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 1) * 
        ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
            1) * log((-log(1 - (1 - (1 - p2))^teta))) * (delta * 
            ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - 
                p2))^teta))) + (-log(1 - (1 - (1 - p2))^teta))^(delta - 
            1) * ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - 
            (1 - p2))^teta))) - 1/delta^2 * ((-log(1 - (1 - (1 - 
            p2))^teta))^(delta - 1) * (delta * ((1 - (1 - p2))^(teta - 
            1) * teta/(1 - (1 - (1 - p2))^teta)))))) - exp(-((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
        (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
            (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
            ((-log(1 - (1 - (1 - p1))^teta))^delta * log((-log(1 - 
                (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - 
                p2))^teta))^delta * log((-log(1 - (1 - (1 - p2))^teta))))) - 
            ((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - 
                (1 - (1 - p2))^teta))^delta)^(1/delta) * (log(((-log(1 - 
                (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
                (1 - p2))^teta))^delta)) * (1/delta^2))) * (((-log(1 - 
        (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
        1) * ((1/delta) * ((-log(1 - (1 - (1 - p2))^teta))^(delta - 
        1) * (delta * ((1 - (1 - p2))^(teta - 1) * teta/(1 - 
        (1 - (1 - p2))^teta)))))))) - exp(-((-log(1 - (1 - (1 - 
    p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (exp(-((-log(1 - (1 - (1 - p1))^teta))^delta + 
    (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * (((-log(1 - 
    (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^((1/delta) - 
    1) * ((1/delta) * ((-log(1 - (1 - (1 - p1))^teta))^delta * 
    log((-log(1 - (1 - (1 - p1))^teta))) + (-log(1 - (1 - (1 - 
    p2))^teta))^delta * log((-log(1 - (1 - (1 - p2))^teta))))) - 
    ((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - (1 - 
        p2))^teta))^delta)^(1/delta) * (log(((-log(1 - (1 - (1 - 
        p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)) * 
        (1/delta^2))))) * ((1/teta) * (exp(-((-log(1 - (1 - (1 - 
    p1))^teta))^delta + (-log(1 - (1 - (1 - p2))^teta))^delta)^(1/delta)) * 
    (((-log(1 - (1 - (1 - p1))^teta))^delta + (-log(1 - (1 - 
        (1 - p2))^teta))^delta)^((1/delta) - 1) * ((1/delta) * 
        ((-log(1 - (1 - (1 - p2))^teta))^(delta - 1) * (delta * 
            ((1 - (1 - p2))^(teta - 1) * teta/(1 - (1 - (1 - 
                p2))^teta))))))))))*exp(delta.st)





bit1.thdel <--(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^((1/(exp(teta.st) + 1)) - 1) * 
    ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        ((((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^(((1/(exp(delta.st) + 1)) - 1) - 1) * (((1/(exp(delta.st) + 
            1)) - 1) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
            p1))^(exp(teta.st) + 1)))) * exp(delta.st)) + (-log(1 - 
            (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))) * exp(delta.st)))) - ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
            1) * (log(((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
            1)^2))) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))))) + ((-log(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
            (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * (log((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) * exp(delta.st)) * ((exp(delta.st) + 1) * ((1 - 
            (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * 
            exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * (exp(delta.st) * 
            ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - 
                p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
                1)))) + ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * (log((-log(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))) * exp(delta.st)) * 
            ((exp(delta.st) + 1) * ((1 - (1 - p2))^(exp(teta.st) + 
                1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - 
                (1 - (1 - p2))^(exp(teta.st) + 1)))) + (-log(1 - 
            (1 - (1 - p2))^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * (exp(delta.st) * ((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))))) - exp(delta.st)/(exp(delta.st) + 
            1)^2 * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
            1) * ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
            (1 - p1))) * exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
            1) * ((1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - 
            (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1))))))) - exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
            p2))^(exp(teta.st) + 1)))) * exp(delta.st)))) - ((-log(1 - 
            (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
            (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * 
        (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (1 - p1))) * exp(teta.st))/(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))) + (-log(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))^((exp(delta.st) + 
            1) - 1) * ((exp(delta.st) + 1) * ((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1))))))))) - exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)))^(((1/(exp(teta.st) + 1)) - 
    1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - 
    p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))) * exp(delta.st)) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) * 
    (log((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))) * exp(delta.st)))) - 
    ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - (1 - (1 - 
        p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))))) * ((1/(exp(teta.st) + 
    1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 1) * 
    ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * 
        ((1 - (1 - p1))^(exp(teta.st) + 1) * (log((1 - (1 - p1))) * 
            exp(teta.st))/(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * 
        ((1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - (1 - p2))) * 
            exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1))))))))) - (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)))^(1/(exp(teta.st) + 
    1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 1) * 
    ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))) * exp(delta.st)))) - ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)))/exp(-((-log(1 - 
    (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^((1/(exp(teta.st) + 1)) - 1) * 
    ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        (((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
            1)) * ((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 
            1)))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
            p2))^(exp(teta.st) + 1)))) * exp(delta.st)))) - ((-log(1 - 
            (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - (1 - 
            (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
            (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))))) * 
    (log(exp(-((-log(1 - (1 - (1 - p1))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2))))



}


if(BivD=="BB6.270"){

  
   
  c.copula.be1 <- 1 - exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
    (1 - p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
        1) * (delta * ((1 - p1)^-(teta + 1) * teta/(1 - (1 - 
        p1)^-teta))))))))



 
  c.copula.be2 <- exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p2))^-teta))))))))

  c.copula.theta <- (-(exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^(-1/teta) * (log(exp(-((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta))) * 
    (1/teta^2)) - exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
    (1 - (1 - p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * 
    ((-1/teta) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
        (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
        1) * (delta * ((1 - p1)^-teta * log((1 - p1))/(1 - (1 - 
        p1)^-teta))) + (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^-teta * log((1 - (1 - p2)))/(1 - 
        (1 - (1 - p2))^-teta))))))))))*(-exp(teta.st))


  c.copula.delta <- (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta) * 
        (log(((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)) * (1/delta^2)) - ((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-delta * 
        log((-log(1 - (1 - (1 - p2))^-teta))) + (-log(1 - (1 - 
        p1)^-teta))^-delta * log((-log(1 - (1 - p1)^-teta)))))))))*(-exp(delta.st))


 

 c.copula2.be1 <-  -(exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^(((-1/teta) - 1) - 1) * 
    (((-1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
        (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * 
                ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta)))))))) * 
    ((-1/teta) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
        (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
        1) * (delta * ((1 - p1)^-(teta + 1) * teta/(1 - (1 - 
        p1)^-teta)))))))) + exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
    1) * (delta * ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta)))))) * 
    (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
        p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * ((1 - 
            p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta)))))) + 
    exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
        1) * (delta * ((1 - p1)^-(teta + 1 + 1) * (teta + 1) * 
        teta/(1 - (1 - p1)^-teta) + (1 - p1)^-(teta + 1) * teta * 
        ((1 - p1)^-(teta + 1) * teta)/(1 - (1 - p1)^-teta)^2)) - 
        (-log(1 - (1 - p1)^-teta))^-(delta + 1 + 1) * ((delta + 
            1) * ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta))) * 
            (delta * ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta))))) - 
        ((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
            p2))^-teta))^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
            1) * ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * 
            ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta))))) * 
            ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
                1) * (delta * ((1 - p1)^-(teta + 1) * teta/(1 - 
                (1 - p1)^-teta)))))))))


                  
 c.copula2.be2 <- exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(((-1/delta) - 
        1) - 1) * (((-1/delta) - 1) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p2))^-teta))))) * ((-1/delta) * ((-log(1 - 
        (1 - (1 - p2))^-teta))^-(delta + 1) * (delta * ((1 - 
        (1 - p2))^-(teta + 1) * teta/(1 - (1 - (1 - p2))^-teta))))) + 
        ((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
            p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1 + 1) * 
                ((delta + 1) * ((1 - (1 - p2))^-(teta + 1) * 
                  teta/(1 - (1 - (1 - p2))^-teta))) * (delta * 
                ((1 - (1 - p2))^-(teta + 1) * teta/(1 - (1 - 
                  (1 - p2))^-teta))) - (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
                1) * (delta * ((1 - (1 - p2))^-(teta + 1 + 1) * 
                (teta + 1) * teta/(1 - (1 - (1 - p2))^-teta) + 
                (1 - (1 - p2))^-(teta + 1) * teta * ((1 - (1 - 
                  p2))^-(teta + 1) * teta)/(1 - (1 - (1 - p2))^-teta)^2))))) - 
        exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
            (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
            1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
            (1 - (1 - p2))^-teta)))))) * (((-log(1 - (1 - p1)^-teta))^-delta + 
            (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
            1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
            (1 - (1 - p2))^-teta)))))))) - exp(-((-log(1 - (1 - 
    p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
    1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - (1 - 
    (1 - p2))^-teta)))))))) * ((-1/teta) * (exp(-((-log(1 - (1 - 
    p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
        p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1) * (delta * 
            ((1 - (1 - p2))^-(teta + 1) * teta/(1 - (1 - (1 - 
                p2))^-teta))))))))

                  


c.copula2.be1be2 <- -(exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(((-1/delta) - 
        1) - 1) * (((-1/delta) - 1) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p2))^-teta))))) * ((-1/delta) * ((-log(1 - 
        (1 - p1)^-teta))^-(delta + 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta/(1 - (1 - p1)^-teta)))))) - exp(-((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1) * 
                (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
                  (1 - (1 - p2))^-teta)))))) * (((-log(1 - (1 - 
        p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
        1) * (delta * ((1 - p1)^-(teta + 1) * teta/(1 - (1 - 
        p1)^-teta)))))))) - exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
    1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - (1 - 
    (1 - p2))^-teta)))))))) * ((-1/teta) * (exp(-((-log(1 - (1 - 
    p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
        p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * ((1 - 
            p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta)))))))))




c.copula2.be1th <-(-((exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * (log(exp(-((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta))) * 
    (1/teta^2)) - exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
    (1 - (1 - p2))^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
    1) * (delta * ((1 - p1)^-teta * log((1 - p1))/(1 - (1 - p1)^-teta))) + 
    (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1) * (delta * 
        ((1 - (1 - p2))^-teta * log((1 - (1 - p2)))/(1 - (1 - 
            (1 - p2))^-teta))))))))) * ((-1/teta) * (exp(-((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
        p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * ((1 - 
            p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta)))))))) + 
    exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * 
        (1/teta^2 * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
            (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
            (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
                (1 - p2))^-teta))^-delta)^((-1/delta) - 1) * 
                ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
                  1) * (delta * ((1 - p1)^-(teta + 1) * teta/(1 - 
                  (1 - p1)^-teta))))))) + (-1/teta) * (exp(-((-log(1 - 
            (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
            (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
                (1 - p2))^-teta))^-delta)^(((-1/delta) - 1) - 
                1) * (((-1/delta) - 1) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
                1) * (delta * ((1 - p1)^-teta * log((1 - p1))/(1 - 
                (1 - p1)^-teta))) + (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
                1) * (delta * ((1 - (1 - p2))^-teta * log((1 - 
                (1 - p2)))/(1 - (1 - (1 - p2))^-teta))))) * ((-1/delta) * 
                ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * 
                  ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta))))) + 
                ((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
                  (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
                  1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
                  1 + 1) * ((delta + 1) * ((1 - p1)^-teta * log((1 - 
                  p1))/(1 - (1 - p1)^-teta))) * (delta * ((1 - 
                  p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta))) + 
                  (-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * 
                    (((1 - p1)^-(teta + 1) - (1 - p1)^-(teta + 
                      1) * log((1 - p1)) * teta)/(1 - (1 - p1)^-teta) - 
                      (1 - p1)^-(teta + 1) * teta * ((1 - p1)^-teta * 
                        log((1 - p1)))/(1 - (1 - p1)^-teta)^2))))) - 
            exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
                (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
                (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
                  (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
                  1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
                  1) * (delta * ((1 - p1)^-teta * log((1 - p1))/(1 - 
                  (1 - p1)^-teta))) + (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
                  1) * (delta * ((1 - (1 - p2))^-teta * log((1 - 
                  (1 - p2)))/(1 - (1 - (1 - p2))^-teta)))))) * 
                (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
                  (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
                  1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
                  1) * (delta * ((1 - p1)^-(teta + 1) * teta/(1 - 
                  (1 - p1)^-teta))))))))))*(-exp(teta.st))




c.copula2.be2th <-((exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * (log(exp(-((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta))) * 
    (1/teta^2)) - exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
    (1 - (1 - p2))^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
    1) * (delta * ((1 - p1)^-teta * log((1 - p1))/(1 - (1 - p1)^-teta))) + 
    (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1) * (delta * 
        ((1 - (1 - p2))^-teta * log((1 - (1 - p2)))/(1 - (1 - 
            (1 - p2))^-teta))))))))) * ((-1/teta) * (exp(-((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
        p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1) * (delta * 
            ((1 - (1 - p2))^-(teta + 1) * teta/(1 - (1 - (1 - 
                p2))^-teta)))))))) + exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 
    1) * (1/teta^2 * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
    1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - (1 - 
    (1 - p2))^-teta))))))) + (-1/teta) * (exp(-((-log(1 - (1 - 
    p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
        p2))^-teta))^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
        1) * ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * 
        ((1 - p1)^-teta * log((1 - p1))/(1 - (1 - p1)^-teta))) + 
        (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1) * (delta * 
            ((1 - (1 - p2))^-teta * log((1 - (1 - p2)))/(1 - 
                (1 - (1 - p2))^-teta))))) * ((-1/delta) * ((-log(1 - 
        (1 - (1 - p2))^-teta))^-(delta + 1) * (delta * ((1 - 
        (1 - p2))^-(teta + 1) * teta/(1 - (1 - (1 - p2))^-teta))))) + 
        ((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
            p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1 + 1) * 
                ((delta + 1) * ((1 - (1 - p2))^-teta * log((1 - 
                  (1 - p2)))/(1 - (1 - (1 - p2))^-teta))) * (delta * 
                ((1 - (1 - p2))^-(teta + 1) * teta/(1 - (1 - 
                  (1 - p2))^-teta))) + (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
                1) * (delta * (((1 - (1 - p2))^-(teta + 1) - 
                (1 - (1 - p2))^-(teta + 1) * log((1 - (1 - p2))) * 
                  teta)/(1 - (1 - (1 - p2))^-teta) - (1 - (1 - 
                p2))^-(teta + 1) * teta * ((1 - (1 - p2))^-teta * 
                log((1 - (1 - p2))))/(1 - (1 - (1 - p2))^-teta)^2))))) - 
    exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - (1 - 
        p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
        1) * (delta * ((1 - p1)^-teta * log((1 - p1))/(1 - (1 - 
        p1)^-teta))) + (-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^-teta * log((1 - (1 - p2)))/(1 - 
        (1 - (1 - p2))^-teta)))))) * (((-log(1 - (1 - p1)^-teta))^-delta + 
        (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p2))^-teta)))))))))*(-exp(teta.st))







bit1.th2 <--((exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
    ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - 
        (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-(delta + 
        1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
        p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 1))))))))) - 
    exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))^((1/(exp(teta.st) + 
        1)) - 1) * (log(exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^-delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^-delta)^(-1/delta))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2))) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
    ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - 
        (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-(delta + 
        1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
        p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 1))))))))) + 
    exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))^((1/(exp(teta.st) + 
        1)) - 1) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - 
        (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - 
        (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - 
            (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^-(delta + 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))) + (-log(1 - (1 - 
            p1)^(exp(teta.st) + 1)))^-(delta + 1) * (delta * 
            ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st))/(1 - 
                (1 - p1)^(exp(teta.st) + 1))))))) * (((-log(1 - 
        (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - 
        (1 - p2))^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1)))) + (-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
        1))))))) + exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + 
        (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - 
            (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^-(delta + 1) * (delta * (((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st)) * (log((1 - 
            (1 - p2))) * exp(teta.st)) + (1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st)))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)) + (1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st)) * ((1 - 
            (1 - p2))^(exp(teta.st) + 1) * (log((1 - (1 - p2))) * 
            exp(teta.st)))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1))^2)) - (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^-(delta + 1 + 1) * ((delta + 1) * ((1 - (1 - 
            p2))^(exp(teta.st) + 1) * (log((1 - (1 - p2))) * 
            exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))) + ((-log(1 - (1 - 
            p1)^(exp(teta.st) + 1)))^-(delta + 1) * (delta * 
            (((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
                exp(teta.st)) * (log((1 - p1)) * exp(teta.st)) + 
                (1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
                  exp(teta.st)))/(1 - (1 - p1)^(exp(teta.st) + 
                1)) + (1 - p1)^(exp(teta.st) + 1) * (log((1 - 
                p1)) * exp(teta.st)) * ((1 - p1)^(exp(teta.st) + 
                1) * (log((1 - p1)) * exp(teta.st)))/(1 - (1 - 
                p1)^(exp(teta.st) + 1))^2)) - (-log(1 - (1 - 
            p1)^(exp(teta.st) + 1)))^-(delta + 1 + 1) * ((delta + 
            1) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
            exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 1)))) * 
            (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
                p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
                1))))))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
            1)))^-delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^-delta)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
            1) * ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-(delta + 
            1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 1) * 
            (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
            p2))^(exp(teta.st) + 1)))) + (-log(1 - (1 - p1)^(exp(teta.st) + 
            1)))^-(delta + 1) * (delta * ((1 - p1)^(exp(teta.st) + 
            1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
            1)))))) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^-(delta + 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))) + (-log(1 - (1 - 
            p1)^(exp(teta.st) + 1)))^-(delta + 1) * (delta * 
            ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st))/(1 - 
                (1 - p1)^(exp(teta.st) + 1)))))))) - exp(teta.st)/(exp(teta.st) + 
        1)^2 * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + 
        (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - 
            (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^-(delta + 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1)))) + (-log(1 - (1 - 
            p1)^(exp(teta.st) + 1)))^-(delta + 1) * (delta * 
            ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st))/(1 - 
                (1 - p1)^(exp(teta.st) + 1))))))))) - ((exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - 
    p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^-delta)^(-1/delta)) * (((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^-delta + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^-(delta + 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
    1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))) + (-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^-(delta + 1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * 
    (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
    1))))))))) - exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))^(1/(exp(teta.st) + 
    1)) * (log(exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2))) * (log(exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))^(1/(exp(teta.st) + 
    1)) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^((-1/delta) - 
        1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1)))) + (-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^-(delta + 1) * (delta * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
        1)))))))/exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^-delta + 
    (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta)) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2) + log(exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^-delta + (-log(1 - (1 - (1 - 
    p2))^(exp(teta.st) + 1)))^-delta)^(-1/delta))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
    1)))/((exp(teta.st) + 1)^2)^2))))





bit1.del2 <-exp(-((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^((-1/teta) - 1) * ((-1/teta) * (exp(-((-log(1 - (1 - 
    p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * ((((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^(((1/(exp(delta.st) + 
    1)) - 1) - 1) * (((1/(exp(delta.st) + 1)) - 1) * ((-log(1 - 
    (1 - p1)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    p1)^-teta))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p2))^-teta))) * exp(delta.st)))) - 
    ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * (log(((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * ((1/(exp(delta.st) + 
    1)) * ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - p1)^-teta))) * exp(delta.st)) + (-log(1 - (1 - (1 - 
    p2))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
    p2))^-teta))) * exp(delta.st)))) + ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^-teta))) * exp(delta.st)) * 
    (log((-log(1 - (1 - p1)^-teta))) * exp(delta.st)) + (-log(1 - 
    (1 - p1)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
    p1)^-teta))) * exp(delta.st)) + ((-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - (1 - p2))^-teta))) * exp(delta.st)) * 
    (log((-log(1 - (1 - (1 - p2))^-teta))) * exp(delta.st)) + 
    (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^-teta))) * exp(delta.st)))) - exp(delta.st)/(exp(delta.st) + 
    1)^2 * ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) * 
    (log((-log(1 - (1 - p1)^-teta))) * exp(delta.st)) + (-log(1 - 
    (1 - (1 - p2))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
    (1 - (1 - p2))^-teta))) * exp(delta.st)))) - ((((-log(1 - 
    (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - 
    p2))^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
    1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^-teta))) * exp(delta.st)) + 
    (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^-teta))) * exp(delta.st)))) - ((-log(1 - 
    (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - 
    p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
    (log(((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2))) * (log(((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))) * 
    (exp(delta.st)/(exp(delta.st) + 1)^2)) + ((-log(1 - (1 - 
    p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)) * (((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^-teta))) * exp(delta.st)) + 
    (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^-teta))) * exp(delta.st)))/((-log(1 - 
    (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - 
    p2))^-teta))^(exp(delta.st) + 1)) * (exp(delta.st)/(exp(delta.st) + 
    1)^2) + log(((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))) * 
    (exp(delta.st)/(exp(delta.st) + 1)^2 - exp(delta.st) * (2 * 
        (exp(delta.st) * (exp(delta.st) + 1)))/((exp(delta.st) + 
        1)^2)^2)))) - exp(-((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) + 
    (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^-teta))) * exp(delta.st)) + 
    (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^-teta))) * exp(delta.st)))) - ((-log(1 - 
    (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - 
    p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
    (log(((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2))) * (((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) * (log((-log(1 - (1 - p1)^-teta))) * exp(delta.st)) + 
    (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^-teta))) * exp(delta.st)))) - ((-log(1 - 
    (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - 
    p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * 
    (log(((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))) * (exp(delta.st)/(exp(delta.st) + 
        1)^2))))) - exp(-((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^(((-1/teta) - 1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - 
    (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - 
    p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
    (((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
        1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - 
        p1)^-teta))^(exp(delta.st) + 1) * (log((-log(1 - (1 - 
        p1)^-teta))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
        1) * (log((-log(1 - (1 - (1 - p2))^-teta))) * exp(delta.st)))) - 
        ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - 
            (1 - (1 - p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
            1)) * (log(((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
            1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
            1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))))) * 
    ((-1/teta) * (exp(-((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * ((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 1) * 
        (log((-log(1 - (1 - p1)^-teta))) * exp(delta.st)) + (-log(1 - 
        (1 - (1 - p2))^-teta))^(exp(delta.st) + 1) * (log((-log(1 - 
        (1 - (1 - p2))^-teta))) * exp(delta.st)))) - ((-log(1 - 
        (1 - p1)^-teta))^(exp(delta.st) + 1) + (-log(1 - (1 - 
        (1 - p2))^-teta))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
        1)) * (log(((-log(1 - (1 - p1)^-teta))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^-teta))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2)))))



c.copula2.be1del <-(-(exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * ((((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * (log(((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
        (1 - (1 - p2))^-teta))^-delta)) * (1/delta^2)) - ((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(((-1/delta) - 
        1) - 1) * (((-1/delta) - 1) * ((-log(1 - (1 - (1 - p2))^-teta))^-delta * 
        log((-log(1 - (1 - (1 - p2))^-teta))) + (-log(1 - (1 - 
        p1)^-teta))^-delta * log((-log(1 - (1 - p1)^-teta)))))) * 
        ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * 
            (delta * ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta))))) + 
        ((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
            p2))^-teta))^-delta)^((-1/delta) - 1) * (1/delta^2 * 
            ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * 
                ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta)))) + 
            (-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
                1) * ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta)) - 
                (-log(1 - (1 - p1)^-teta))^-(delta + 1) * log((-log(1 - 
                  (1 - p1)^-teta))) * (delta * ((1 - p1)^-(teta + 
                  1) * teta/(1 - (1 - p1)^-teta)))))) - exp(-((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)^(-1/delta) * (log(((-log(1 - 
            (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)) * 
            (1/delta^2)) - ((-log(1 - (1 - p1)^-teta))^-delta + 
            (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-delta * 
            log((-log(1 - (1 - (1 - p2))^-teta))) + (-log(1 - 
            (1 - p1)^-teta))^-delta * log((-log(1 - (1 - p1)^-teta)))))) * 
        (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - p1)^-teta))^-(delta + 1) * (delta * 
                ((1 - p1)^-(teta + 1) * teta/(1 - (1 - p1)^-teta)))))))) - 
    exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta))^(((-1/teta) - 1) - 
        1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
        (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)^(-1/delta) * (log(((-log(1 - 
            (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)) * 
            (1/delta^2)) - ((-log(1 - (1 - p1)^-teta))^-delta + 
            (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-delta * 
            log((-log(1 - (1 - (1 - p2))^-teta))) + (-log(1 - 
            (1 - p1)^-teta))^-delta * log((-log(1 - (1 - p1)^-teta)))))))) * 
        ((-1/teta) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
            (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
            (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
                (1 - p2))^-teta))^-delta)^((-1/delta) - 1) * 
                ((-1/delta) * ((-log(1 - (1 - p1)^-teta))^-(delta + 
                  1) * (delta * ((1 - p1)^-(teta + 1) * teta/(1 - 
                  (1 - p1)^-teta))))))))))*(-exp(delta.st))


        



        
c.copula2.be2del <- (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
    p2))^-teta))^-delta)^(-1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (exp(-((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)^(-1/delta)) * ((((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * (log(((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - 
        (1 - (1 - p2))^-teta))^-delta)) * (1/delta^2)) - ((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(((-1/delta) - 
        1) - 1) * (((-1/delta) - 1) * ((-log(1 - (1 - (1 - p2))^-teta))^-delta * 
        log((-log(1 - (1 - (1 - p2))^-teta))) + (-log(1 - (1 - 
        p1)^-teta))^-delta * log((-log(1 - (1 - p1)^-teta)))))) * 
        ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
            1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
            (1 - (1 - p2))^-teta))))) + ((-log(1 - (1 - p1)^-teta))^-delta + 
        (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
        1) * (1/delta^2 * ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 
        1) * (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
        (1 - (1 - p2))^-teta)))) + (-1/delta) * ((-log(1 - (1 - 
        (1 - p2))^-teta))^-(delta + 1) * ((1 - (1 - p2))^-(teta + 
        1) * teta/(1 - (1 - (1 - p2))^-teta)) - (-log(1 - (1 - 
        (1 - p2))^-teta))^-(delta + 1) * log((-log(1 - (1 - (1 - 
        p2))^-teta))) * (delta * ((1 - (1 - p2))^-(teta + 1) * 
        teta/(1 - (1 - (1 - p2))^-teta)))))) - exp(-((-log(1 - 
        (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
        (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)^(-1/delta) * (log(((-log(1 - 
            (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)) * 
            (1/delta^2)) - ((-log(1 - (1 - p1)^-teta))^-delta + 
            (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
            1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-delta * 
            log((-log(1 - (1 - (1 - p2))^-teta))) + (-log(1 - 
            (1 - p1)^-teta))^-delta * log((-log(1 - (1 - p1)^-teta)))))) * 
        (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
            (1 - p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
            ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1) * 
                (delta * ((1 - (1 - p2))^-(teta + 1) * teta/(1 - 
                  (1 - (1 - p2))^-teta)))))))) - exp(-((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (exp(-((-log(1 - (1 - p1)^-teta))^-delta + 
    (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * (((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta) * 
    (log(((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - 
        (1 - p2))^-teta))^-delta)) * (1/delta^2)) - ((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^((-1/delta) - 
    1) * ((-1/delta) * ((-log(1 - (1 - (1 - p2))^-teta))^-delta * 
    log((-log(1 - (1 - (1 - p2))^-teta))) + (-log(1 - (1 - p1)^-teta))^-delta * 
    log((-log(1 - (1 - p1)^-teta)))))))) * ((-1/teta) * (exp(-((-log(1 - 
    (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - p2))^-teta))^-delta)^(-1/delta)) * 
    (((-log(1 - (1 - p1)^-teta))^-delta + (-log(1 - (1 - (1 - 
        p2))^-teta))^-delta)^((-1/delta) - 1) * ((-1/delta) * 
        ((-log(1 - (1 - (1 - p2))^-teta))^-(delta + 1) * (delta * 
            ((1 - (1 - p2))^-(teta + 1) * teta/(1 - (1 - (1 - 
                p2))^-teta)))))))))*(-exp(delta.st))



bit1.thdel <-exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1)))^((1/(exp(teta.st) + 1)) - 1) * 
    ((1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * 
        ((((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1))^(((1/(exp(delta.st) + 1)) - 1) - 1) * (((1/(exp(delta.st) + 
            1)) - 1) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
            1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 1)))) * 
            exp(delta.st)) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - 
            p2))^(exp(teta.st) + 1)))) * exp(delta.st)))) - ((-log(1 - 
            (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
            (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))^((1/(exp(delta.st) + 1)) - 1) * (log(((-log(1 - 
            (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + 
            (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
                1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * 
            ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
                1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
                1) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
                p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
                1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
                1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 
                1) * ((1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - 
                (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
                1)))))) + ((-log(1 - (1 - p1)^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
            1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * (log((-log(1 - (1 - 
            p1)^(exp(teta.st) + 1)))) * exp(delta.st)) * ((exp(delta.st) + 
            1) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
            exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 1)))) + 
            (-log(1 - (1 - p1)^(exp(teta.st) + 1)))^((exp(delta.st) + 
                1) - 1) * (exp(delta.st) * ((1 - p1)^(exp(teta.st) + 
                1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - 
                p1)^(exp(teta.st) + 1)))) + ((-log(1 - (1 - (1 - 
            p2))^(exp(teta.st) + 1)))^((exp(delta.st) + 1) - 
            1) * (log((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))) * exp(delta.st)) * ((exp(delta.st) + 1) * ((1 - 
            (1 - p2))^(exp(teta.st) + 1) * (log((1 - (1 - p2))) * 
            exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
            1)))^((exp(delta.st) + 1) - 1) * (exp(delta.st) * 
            ((1 - (1 - p2))^(exp(teta.st) + 1) * (log((1 - (1 - 
                p2))) * exp(teta.st))/(1 - (1 - (1 - p2))^(exp(teta.st) + 
                1)))))) - exp(delta.st)/(exp(delta.st) + 1)^2 * 
            ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^((exp(delta.st) + 
                1) - 1) * ((exp(delta.st) + 1) * ((1 - p1)^(exp(teta.st) + 
                1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - 
                p1)^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - 
                p2))^(exp(teta.st) + 1)))^((exp(delta.st) + 1) - 
                1) * ((exp(delta.st) + 1) * ((1 - (1 - p2))^(exp(teta.st) + 
                1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - 
                (1 - (1 - p2))^(exp(teta.st) + 1))))))) - exp(-((-log(1 - 
        (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
        1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - 
        (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))) * (((-log(1 - 
        (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^((1/(exp(delta.st) + 1)) - 1) * ((1/(exp(delta.st) + 
        1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^((exp(delta.st) + 
        1) - 1) * ((exp(delta.st) + 1) * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st))/(1 - (1 - p1)^(exp(teta.st) + 
        1)))) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^((exp(delta.st) + 
        1) - 1) * ((exp(delta.st) + 1) * ((1 - (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1))))))))) - exp(-((-log(1 - (1 - 
    p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^(((1/(exp(teta.st) + 1)) - 1) - 1) * (((1/(exp(teta.st) + 
    1)) - 1) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 1) * 
    ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))) * 
    (exp(delta.st)/(exp(delta.st) + 1)^2))))) * ((1/(exp(teta.st) + 
    1)) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
    1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 1) * 
    ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^((exp(delta.st) + 1) - 1) * ((exp(delta.st) + 1) * 
        ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st))/(1 - 
            (1 - p1)^(exp(teta.st) + 1)))) + (-log(1 - (1 - (1 - 
        p2))^(exp(teta.st) + 1)))^((exp(delta.st) + 1) - 1) * 
        ((exp(delta.st) + 1) * ((1 - (1 - p2))^(exp(teta.st) + 
            1) * (log((1 - (1 - p2))) * exp(teta.st))/(1 - (1 - 
            (1 - p2))^(exp(teta.st) + 1))))))))) - (exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^(1/(exp(teta.st) + 1)) * (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1))) * (((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 
    1)) - 1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
    1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))) * 
    (exp(delta.st)/(exp(delta.st) + 1)^2)))/exp(-((-log(1 - (1 - 
    p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1))) * (exp(teta.st)/(exp(teta.st) + 1)^2)) + exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
    (exp(-((-log(1 - (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))^(1/(exp(delta.st) + 1))) * (((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^((1/(exp(delta.st) + 1)) - 
        1) * ((1/(exp(delta.st) + 1)) * ((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))) * exp(delta.st)) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) * (log((-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))) * exp(delta.st)))) - ((-log(1 - (1 - p1)^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1) + (-log(1 - (1 - (1 - p2))^(exp(teta.st) + 
        1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 1)) * (log(((-log(1 - 
        (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
        (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 
        1))) * (exp(delta.st)/(exp(delta.st) + 1)^2))))) * (log(exp(-((-log(1 - 
    (1 - p1)^(exp(teta.st) + 1)))^(exp(delta.st) + 1) + (-log(1 - 
    (1 - (1 - p2))^(exp(teta.st) + 1)))^(exp(delta.st) + 1))^(1/(exp(delta.st) + 
    1)))) * (exp(teta.st)/(exp(teta.st) + 1)^2)))



}


if(BivD=="BB7.0"){
  
   
  c.copula.be1 <- -((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * (((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
    ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
        1) * teta)))))))

 
  c.copula.be2 <- -((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * (((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
    ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * ((1 - p2)^(teta - 
        1) * teta)))))))

  c.copula.theta <- ((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^(1/teta) * (log((1 - ((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^(-1/delta))) * (1/teta^2)) + 
    (1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * (((1 - (1 - 
        p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 1) * 
        (delta * ((1 - p1)^teta * log((1 - p1)))) + (1 - (1 - 
        p2)^teta)^-(delta + 1) * (delta * ((1 - p2)^teta * log((1 - 
        p2)))))))))*exp(teta.st)



  c.copula.delta <- ((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * (((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^(-1/delta) * (log(((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)) * (1/delta^2)) - 
    ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-delta * 
        log((1 - (1 - p2)^teta)) + (1 - (1 - p1)^teta)^-delta * 
        log((1 - (1 - p1)^teta)))))))*exp(delta.st)




 c.copula2.be1 <- -((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * 
    ((1 - p1)^(teta - 1) * teta)))))) * ((1/teta) * (((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * 
    ((1 - p1)^(teta - 1) * teta)))))) - (1 - ((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^(-1/delta))^((1/teta) - 1) * 
    ((1/teta) * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 
        1) * (delta * ((1 - p1)^((teta - 1) - 1) * (teta - 1) * 
        teta)) + (1 - (1 - p1)^teta)^-(delta + 1 + 1) * ((delta + 
        1) * ((1 - p1)^(teta - 1) * teta)) * (delta * ((1 - p1)^(teta - 
        1) * teta)))) + ((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
        p2)^teta)^-delta - 1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
        1) * ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * ((1 - 
        p1)^(teta - 1) * teta)))) * ((-1/delta) * ((1 - (1 - 
        p1)^teta)^-(delta + 1) * (delta * ((1 - p1)^(teta - 1) * 
        teta)))))))


                  
 c.copula2.be2 <- -((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * 
    ((1 - p2)^(teta - 1) * teta)))))) * ((1/teta) * (((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * 
    ((1 - p2)^(teta - 1) * teta)))))) - (1 - ((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^(-1/delta))^((1/teta) - 1) * 
    ((1/teta) * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-(delta + 
        1) * (delta * ((1 - p2)^((teta - 1) - 1) * (teta - 1) * 
        teta)) + (1 - (1 - p2)^teta)^-(delta + 1 + 1) * ((delta + 
        1) * ((1 - p2)^(teta - 1) * teta)) * (delta * ((1 - p2)^(teta - 
        1) * teta)))) + ((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
        p2)^teta)^-delta - 1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
        1) * ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * ((1 - 
        p2)^(teta - 1) * teta)))) * ((-1/delta) * ((1 - (1 - 
        p2)^teta)^-(delta + 1) * (delta * ((1 - p2)^(teta - 1) * 
        teta)))))))



                 


c.copula2.be1be2 <--((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * 
    ((1 - p2)^(teta - 1) * teta)))))) * ((1/teta) * (((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * 
    ((1 - p1)^(teta - 1) * teta)))))) - (1 - ((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^(-1/delta))^((1/teta) - 1) * 
    ((1/teta) * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((1 - 
        (1 - p2)^teta)^-(delta + 1) * (delta * ((1 - p2)^(teta - 
        1) * teta)))) * ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 
        1) * (delta * ((1 - p1)^(teta - 1) * teta)))))))





c.copula2.be1th <-(-((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * (((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^(((-1/delta) - 1) - 1) * 
    (((-1/delta) - 1) * ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * 
        ((1 - p1)^teta * log((1 - p1)))) + (1 - (1 - p2)^teta)^-(delta + 
        1) * (delta * ((1 - p2)^teta * log((1 - p2)))))) * ((-1/delta) * 
    ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
        1) * teta)))) + ((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
    p2)^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * ((1 - 
    (1 - p1)^teta)^-(delta + 1 + 1) * ((delta + 1) * ((1 - p1)^teta * 
    log((1 - p1)))) * (delta * ((1 - p1)^(teta - 1) * teta)) + 
    (1 - (1 - p1)^teta)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
        1) * log((1 - p1)) * teta + (1 - p1)^(teta - 1)))))) - 
    1/teta^2 * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 
        1) * (delta * ((1 - p1)^(teta - 1) * teta)))))) - ((1 - 
    ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^(-1/delta))^((1/teta) - 1) * (log((1 - ((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^(-1/delta))) * 
    (1/teta^2)) + (1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
    p2)^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 1) - 1) * 
    (((1/teta) - 1) * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
        p2)^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
        ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * ((1 - p1)^teta * 
            log((1 - p1)))) + (1 - (1 - p2)^teta)^-(delta + 1) * 
            (delta * ((1 - p2)^teta * log((1 - p2))))))))) * 
    ((1/teta) * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 
        1) * (delta * ((1 - p1)^(teta - 1) * teta))))))))*exp(teta.st)




c.copula2.be2th <-(-((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * (((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^(((-1/delta) - 1) - 1) * 
    (((-1/delta) - 1) * ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * 
        ((1 - p1)^teta * log((1 - p1)))) + (1 - (1 - p2)^teta)^-(delta + 
        1) * (delta * ((1 - p2)^teta * log((1 - p2)))))) * ((-1/delta) * 
    ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * ((1 - p2)^(teta - 
        1) * teta)))) + ((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
    p2)^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * ((1 - 
    (1 - p2)^teta)^-(delta + 1 + 1) * ((delta + 1) * ((1 - p2)^teta * 
    log((1 - p2)))) * (delta * ((1 - p2)^(teta - 1) * teta)) + 
    (1 - (1 - p2)^teta)^-(delta + 1) * (delta * ((1 - p2)^(teta - 
        1) * log((1 - p2)) * teta + (1 - p2)^(teta - 1)))))) - 
    1/teta^2 * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-(delta + 
        1) * (delta * ((1 - p2)^(teta - 1) * teta)))))) - ((1 - 
    ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^(-1/delta))^((1/teta) - 1) * (log((1 - ((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^(-1/delta))) * 
    (1/teta^2)) + (1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
    p2)^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 1) - 1) * 
    (((1/teta) - 1) * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
        p2)^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
        ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * ((1 - p1)^teta * 
            log((1 - p1)))) + (1 - (1 - p2)^teta)^-(delta + 1) * 
            (delta * ((1 - p2)^teta * log((1 - p2))))))))) * 
    ((1/teta) * (((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-(delta + 
        1) * (delta * ((1 - p2)^(teta - 1) * teta))))))))*exp(teta.st)







bit1.th2 <-(1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-delta - 1)^(-1/delta))^(1/(exp(teta.st) + 1)) * (log((1 - 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 
        1))^-delta - 1)^(-1/delta))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
    1)))/((exp(teta.st) + 1)^2)^2) - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p1)^(exp(teta.st) + 1))^-(delta + 
    1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
    exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(delta + 
    1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
    exp(teta.st))))))/(1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-delta + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^(-1/delta)) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) - ((1 - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))^(1/(exp(teta.st) + 
    1)) * (log((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-delta + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) + (1 - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - p1)^(exp(teta.st) + 
    1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p1)^(exp(teta.st) + 1))^-(delta + 
    1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
    exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(delta + 
    1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
    exp(teta.st))))))))) * (log((1 - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) + ((1 - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - p1)^(exp(teta.st) + 
    1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^(((-1/delta) - 
    1) - 1) * (((-1/delta) - 1) * ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-(delta + 1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * 
    (log((1 - p1)) * exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-(delta + 1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * 
    (log((1 - p2)) * exp(teta.st)))))) * ((-1/delta) * ((1 - 
    (1 - p1)^(exp(teta.st) + 1))^-(delta + 1) * (delta * ((1 - 
    p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st)))) + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^-(delta + 1) * (delta * 
        ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)))))) + 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 
        1))^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * ((1 - 
        (1 - p1)^(exp(teta.st) + 1))^-(delta + 1 + 1) * ((delta + 
        1) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
        exp(teta.st)))) * (delta * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st)))) + (1 - (1 - p1)^(exp(teta.st) + 
        1))^-(delta + 1) * (delta * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st)) * (log((1 - p1)) * 
        exp(teta.st)) + (1 - p1)^(exp(teta.st) + 1) * (log((1 - 
        p1)) * exp(teta.st)))) + ((1 - (1 - p2)^(exp(teta.st) + 
        1))^-(delta + 1 + 1) * ((delta + 1) * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st)))) * (delta * ((1 - 
        p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)))) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-(delta + 1) * (delta * 
            ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)) * 
                (log((1 - p2)) * exp(teta.st)) + (1 - p2)^(exp(teta.st) + 
                1) * (log((1 - p2)) * exp(teta.st)))))))) - exp(teta.st)/(exp(teta.st) + 
    1)^2 * (((1 - (1 - p1)^(exp(teta.st) + 1))^-delta + (1 - 
    (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^((-1/delta) - 1) * 
    ((-1/delta) * ((1 - (1 - p1)^(exp(teta.st) + 1))^-(delta + 
        1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
        p1)) * exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 
        1))^-(delta + 1) * (delta * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st)))))))) - ((1 - ((1 - 
    (1 - p1)^(exp(teta.st) + 1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-delta - 1)^(-1/delta))^((1/(exp(teta.st) + 1)) - 1) * 
    (log((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-delta + (1 - 
        (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))) * 
        (exp(teta.st)/(exp(teta.st) + 1)^2)) + (1 - ((1 - (1 - 
    p1)^(exp(teta.st) + 1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-delta - 1)^(-1/delta))^(((1/(exp(teta.st) + 1)) - 1) - 
    1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 - p1)^(exp(teta.st) + 
    1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 1))^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p1)^(exp(teta.st) + 1))^-(delta + 
    1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
    exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(delta + 
    1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
    exp(teta.st))))))))) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - 
    p1)^(exp(teta.st) + 1))^-delta + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * ((1 - (1 - 
    p1)^(exp(teta.st) + 1))^-(delta + 1) * (delta * ((1 - p1)^(exp(teta.st) + 
    1) * (log((1 - p1)) * exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-(delta + 1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * 
    (log((1 - p2)) * exp(teta.st)))))))))




bit1.del2 <- (1 - ((1 - (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
    1)^(-1/exp(delta.st)))^((1/teta) - 1) * ((1/teta) * ((((1 - 
    (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
    1)^(-1/exp(delta.st)) * (log(((1 - (1 - p1)^teta)^-exp(delta.st) + 
    (1 - (1 - p2)^teta)^-exp(delta.st) - 1)) * (exp(delta.st)/exp(delta.st)^2)) - 
    ((1 - (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
        1)^((-1/exp(delta.st)) - 1) * ((-1/exp(delta.st)) * ((1 - 
        (1 - p2)^teta)^-exp(delta.st) * (log((1 - (1 - p2)^teta)) * 
        exp(delta.st)) + (1 - (1 - p1)^teta)^-exp(delta.st) * 
        (log((1 - (1 - p1)^teta)) * exp(delta.st))))) * (log(((1 - 
    (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
    1)) * (exp(delta.st)/exp(delta.st)^2)) + ((1 - (1 - p1)^teta)^-exp(delta.st) + 
    (1 - (1 - p2)^teta)^-exp(delta.st) - 1)^(-1/exp(delta.st)) * 
    (log(((1 - (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
        1)) * (exp(delta.st)/exp(delta.st)^2 - exp(delta.st) * 
        (2 * (exp(delta.st) * exp(delta.st)))/(exp(delta.st)^2)^2) - 
        ((1 - (1 - p2)^teta)^-exp(delta.st) * (log((1 - (1 - 
            p2)^teta)) * exp(delta.st)) + (1 - (1 - p1)^teta)^-exp(delta.st) * 
            (log((1 - (1 - p1)^teta)) * exp(delta.st)))/((1 - 
            (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
            1) * (exp(delta.st)/exp(delta.st)^2)) - ((((1 - (1 - 
    p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
    1)^((-1/exp(delta.st)) - 1) * (log(((1 - (1 - p1)^teta)^-exp(delta.st) + 
    (1 - (1 - p2)^teta)^-exp(delta.st) - 1)) * (exp(delta.st)/exp(delta.st)^2)) - 
    ((1 - (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
        1)^(((-1/exp(delta.st)) - 1) - 1) * (((-1/exp(delta.st)) - 
        1) * ((1 - (1 - p2)^teta)^-exp(delta.st) * (log((1 - 
        (1 - p2)^teta)) * exp(delta.st)) + (1 - (1 - p1)^teta)^-exp(delta.st) * 
        (log((1 - (1 - p1)^teta)) * exp(delta.st))))) * ((-1/exp(delta.st)) * 
    ((1 - (1 - p2)^teta)^-exp(delta.st) * (log((1 - (1 - p2)^teta)) * 
        exp(delta.st)) + (1 - (1 - p1)^teta)^-exp(delta.st) * 
        (log((1 - (1 - p1)^teta)) * exp(delta.st)))) + ((1 - 
    (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
    1)^((-1/exp(delta.st)) - 1) * (exp(delta.st)/exp(delta.st)^2 * 
    ((1 - (1 - p2)^teta)^-exp(delta.st) * (log((1 - (1 - p2)^teta)) * 
        exp(delta.st)) + (1 - (1 - p1)^teta)^-exp(delta.st) * 
        (log((1 - (1 - p1)^teta)) * exp(delta.st))) + (-1/exp(delta.st)) * 
    ((1 - (1 - p2)^teta)^-exp(delta.st) * (log((1 - (1 - p2)^teta)) * 
        exp(delta.st)) - (1 - (1 - p2)^teta)^-exp(delta.st) * 
        (log((1 - (1 - p2)^teta)) * exp(delta.st)) * (log((1 - 
        (1 - p2)^teta)) * exp(delta.st)) + ((1 - (1 - p1)^teta)^-exp(delta.st) * 
        (log((1 - (1 - p1)^teta)) * exp(delta.st)) - (1 - (1 - 
        p1)^teta)^-exp(delta.st) * (log((1 - (1 - p1)^teta)) * 
        exp(delta.st)) * (log((1 - (1 - p1)^teta)) * exp(delta.st)))))))) - 
    (1 - ((1 - (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
        1)^(-1/exp(delta.st)))^(((1/teta) - 1) - 1) * (((1/teta) - 
        1) * (((1 - (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - 
        p2)^teta)^-exp(delta.st) - 1)^(-1/exp(delta.st)) * (log(((1 - 
        (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
        1)) * (exp(delta.st)/exp(delta.st)^2)) - ((1 - (1 - p1)^teta)^-exp(delta.st) + 
        (1 - (1 - p2)^teta)^-exp(delta.st) - 1)^((-1/exp(delta.st)) - 
        1) * ((-1/exp(delta.st)) * ((1 - (1 - p2)^teta)^-exp(delta.st) * 
        (log((1 - (1 - p2)^teta)) * exp(delta.st)) + (1 - (1 - 
        p1)^teta)^-exp(delta.st) * (log((1 - (1 - p1)^teta)) * 
        exp(delta.st)))))) * ((1/teta) * (((1 - (1 - p1)^teta)^-exp(delta.st) + 
        (1 - (1 - p2)^teta)^-exp(delta.st) - 1)^(-1/exp(delta.st)) * 
        (log(((1 - (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - 
            p2)^teta)^-exp(delta.st) - 1)) * (exp(delta.st)/exp(delta.st)^2)) - 
        ((1 - (1 - p1)^teta)^-exp(delta.st) + (1 - (1 - p2)^teta)^-exp(delta.st) - 
            1)^((-1/exp(delta.st)) - 1) * ((-1/exp(delta.st)) * 
            ((1 - (1 - p2)^teta)^-exp(delta.st) * (log((1 - (1 - 
                p2)^teta)) * exp(delta.st)) + (1 - (1 - p1)^teta)^-exp(delta.st) * 
                (log((1 - (1 - p1)^teta)) * exp(delta.st))))))



c.copula2.be1del <-(-((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * ((((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * (log(((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)) * (1/delta^2)) - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
    p2)^teta)^-delta - 1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
    1) * ((1 - (1 - p2)^teta)^-delta * log((1 - (1 - p2)^teta)) + 
    (1 - (1 - p1)^teta)^-delta * log((1 - (1 - p1)^teta))))) * 
    ((-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * 
        ((1 - p1)^(teta - 1) * teta)))) + ((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 1) * (1/delta^2 * 
    ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
        1) * teta))) + (-1/delta) * ((1 - (1 - p1)^teta)^-(delta + 
    1) * ((1 - p1)^(teta - 1) * teta) - (1 - (1 - p1)^teta)^-(delta + 
    1) * log((1 - (1 - p1)^teta)) * (delta * ((1 - p1)^(teta - 
    1) * teta)))))) - (1 - ((1 - (1 - p1)^teta)^-delta + (1 - 
    (1 - p2)^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 1) - 
    1) * (((1/teta) - 1) * (((1 - (1 - p1)^teta)^-delta + (1 - 
    (1 - p2)^teta)^-delta - 1)^(-1/delta) * (log(((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)) * (1/delta^2)) - ((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-delta * log((1 - 
    (1 - p2)^teta)) + (1 - (1 - p1)^teta)^-delta * log((1 - (1 - 
    p1)^teta)))))) * ((1/teta) * (((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
    ((1 - (1 - p1)^teta)^-(delta + 1) * (delta * ((1 - p1)^(teta - 
        1) * teta))))))))*exp(delta.st)


        


        
c.copula2.be2del <- (-((1 - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * ((((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * (log(((1 - (1 - p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 
    1)) * (1/delta^2)) - ((1 - (1 - p1)^teta)^-delta + (1 - (1 - 
    p2)^teta)^-delta - 1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 
    1) * ((1 - (1 - p2)^teta)^-delta * log((1 - (1 - p2)^teta)) + 
    (1 - (1 - p1)^teta)^-delta * log((1 - (1 - p1)^teta))))) * 
    ((-1/delta) * ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * 
        ((1 - p2)^(teta - 1) * teta)))) + ((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 1) * (1/delta^2 * 
    ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * ((1 - p2)^(teta - 
        1) * teta))) + (-1/delta) * ((1 - (1 - p2)^teta)^-(delta + 
    1) * ((1 - p2)^(teta - 1) * teta) - (1 - (1 - p2)^teta)^-(delta + 
    1) * log((1 - (1 - p2)^teta)) * (delta * ((1 - p2)^(teta - 
    1) * teta)))))) - (1 - ((1 - (1 - p1)^teta)^-delta + (1 - 
    (1 - p2)^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 1) - 
    1) * (((1/teta) - 1) * (((1 - (1 - p1)^teta)^-delta + (1 - 
    (1 - p2)^teta)^-delta - 1)^(-1/delta) * (log(((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)) * (1/delta^2)) - ((1 - (1 - 
    p1)^teta)^-delta + (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - (1 - p2)^teta)^-delta * log((1 - 
    (1 - p2)^teta)) + (1 - (1 - p1)^teta)^-delta * log((1 - (1 - 
    p1)^teta)))))) * ((1/teta) * (((1 - (1 - p1)^teta)^-delta + 
    (1 - (1 - p2)^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
    ((1 - (1 - p2)^teta)^-(delta + 1) * (delta * ((1 - p2)^(teta - 
        1) * teta))))))))*exp(delta.st)



bit1.thdel <-(1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + (1 - 
    (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 1)^(-1/exp(delta.st)))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * ((((1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 
    1)^((-1/exp(delta.st)) - 1) * (log(((1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 
    1)) * (exp(delta.st)/exp(delta.st)^2)) - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 
    1)^(((-1/exp(delta.st)) - 1) - 1) * (((-1/exp(delta.st)) - 
    1) * ((1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) * 
    (log((1 - (1 - p2)^(exp(teta.st) + 1))) * exp(delta.st)) + 
    (1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) * (log((1 - 
        (1 - p1)^(exp(teta.st) + 1))) * exp(delta.st))))) * ((-1/exp(delta.st)) * 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 1) * 
        (exp(delta.st) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
            p1)) * exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 
        1))^-(exp(delta.st) + 1) * (exp(delta.st) * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st)))))) + ((1 - (1 - 
    p1)^(exp(teta.st) + 1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-exp(delta.st) - 1)^((-1/exp(delta.st)) - 1) * (exp(delta.st)/exp(delta.st)^2 * 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 1) * 
        (exp(delta.st) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
            p1)) * exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 
        1))^-(exp(delta.st) + 1) * (exp(delta.st) * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st))))) + (-1/exp(delta.st)) * 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 1) * 
        (exp(delta.st) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
            p1)) * exp(teta.st)))) - (1 - (1 - p1)^(exp(teta.st) + 
        1))^-(exp(delta.st) + 1) * (log((1 - (1 - p1)^(exp(teta.st) + 
        1))) * exp(delta.st)) * (exp(delta.st) * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st)))) + ((1 - (1 - p2)^(exp(teta.st) + 
        1))^-(exp(delta.st) + 1) * (exp(delta.st) * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st)))) - (1 - (1 - p2)^(exp(teta.st) + 
        1))^-(exp(delta.st) + 1) * (log((1 - (1 - p2)^(exp(teta.st) + 
        1))) * exp(delta.st)) * (exp(delta.st) * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st))))))))) - (1 - ((1 - 
    (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-exp(delta.st) - 1)^(-1/exp(delta.st)))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 - 
    p1)^(exp(teta.st) + 1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-exp(delta.st) - 1)^(-1/exp(delta.st)) * (log(((1 - (1 - 
    p1)^(exp(teta.st) + 1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-exp(delta.st) - 1)) * (exp(delta.st)/exp(delta.st)^2)) - 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + (1 - 
        (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 1)^((-1/exp(delta.st)) - 
        1) * ((-1/exp(delta.st)) * ((1 - (1 - p2)^(exp(teta.st) + 
        1))^-exp(delta.st) * (log((1 - (1 - p2)^(exp(teta.st) + 
        1))) * exp(delta.st)) + (1 - (1 - p1)^(exp(teta.st) + 
        1))^-exp(delta.st) * (log((1 - (1 - p1)^(exp(teta.st) + 
        1))) * exp(delta.st)))))) * ((1/(exp(teta.st) + 1)) * 
    (((1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + (1 - 
        (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 1)^((-1/exp(delta.st)) - 
        1) * ((-1/exp(delta.st)) * ((1 - (1 - p1)^(exp(teta.st) + 
        1))^-(exp(delta.st) + 1) * (exp(delta.st) * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 
        1))^-(exp(delta.st) + 1) * (exp(delta.st) * ((1 - p2)^(exp(teta.st) + 
        1) * (log((1 - p2)) * exp(teta.st)))))))) - ((1 - ((1 - 
    (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-exp(delta.st) - 1)^(-1/exp(delta.st)))^(1/(exp(teta.st) + 
    1)) * ((((1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 1)^(-1/exp(delta.st)) * 
    (log(((1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 1)) * 
        (exp(delta.st)/exp(delta.st)^2)) - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 
    1)^((-1/exp(delta.st)) - 1) * ((-1/exp(delta.st)) * ((1 - 
    (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) * (log((1 - (1 - 
    p2)^(exp(teta.st) + 1))) * exp(delta.st)) + (1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) * (log((1 - (1 - p1)^(exp(teta.st) + 1))) * 
    exp(delta.st)))))/(1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 1)^(-1/exp(delta.st))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) + (1 - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 
    1)^(-1/exp(delta.st)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 
    1)) * (((1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 1)^(-1/exp(delta.st)) * 
    (log(((1 - (1 - p1)^(exp(teta.st) + 1))^-exp(delta.st) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 1)) * 
        (exp(delta.st)/exp(delta.st)^2)) - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 
    1)^((-1/exp(delta.st)) - 1) * ((-1/exp(delta.st)) * ((1 - 
    (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) * (log((1 - (1 - 
    p2)^(exp(teta.st) + 1))) * exp(delta.st)) + (1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) * (log((1 - (1 - p1)^(exp(teta.st) + 1))) * 
    exp(delta.st)))))) * (log((1 - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-exp(delta.st) + (1 - (1 - p2)^(exp(teta.st) + 1))^-exp(delta.st) - 
    1)^(-1/exp(delta.st)))) * (exp(teta.st)/(exp(teta.st) + 1)^2)))



}


if(BivD=="BB7.90"){
  
   
  c.copula.be1 <--((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - p1^-teta)^(delta - 
    1) * (delta * (p1^-(teta + 1) * teta)))))))


 
  c.copula.be2 <-1 + (1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^(1/delta))^((-1/teta) - 1) * ((-1/teta) * (((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta))))))

  c.copula.theta <- ((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^(-1/teta) * 
    (log((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
        1)^(1/delta))) * (1/teta^2)) - (1 - ((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 1) * 
    ((-1/teta) * (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
        1)^((1/delta) - 1) * ((1/delta) * ((1 - p1^-teta)^(delta - 
        1) * (delta * (p1^-teta * log(p1))) + (1 - (1 - p2)^-teta)^(delta - 
        1) * (delta * ((1 - p2)^-teta * log((1 - p2)))))))))*(-exp(teta.st))


  c.copula.delta <- (-((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - p1^-teta)^delta * 
    log((1 - p1^-teta)) + (1 - (1 - p2)^-teta)^delta * log((1 - 
    (1 - p2)^-teta)))) - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^(1/delta) * (log(((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)) * (1/delta^2))))))*(-exp(delta.st))



 c.copula2.be1 <- -((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((1 - p1^-teta)^(delta - 
    1) * (delta * (p1^-(teta + 1) * teta)))) * ((1/delta) * ((1 - 
    p1^-teta)^(delta - 1) * (delta * (p1^-(teta + 1) * teta)))) + 
    ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - p1^-teta)^((delta - 1) - 1) * 
        ((delta - 1) * (p1^-(teta + 1) * teta)) * (delta * (p1^-(teta + 
        1) * teta)) - (1 - p1^-teta)^(delta - 1) * (delta * (p1^-(teta + 
        1 + 1) * (teta + 1) * teta)))))) - (1 - ((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((1 - p1^-teta)^delta + (1 - 
    (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - p1^-teta)^(delta - 1) * (delta * (p1^-(teta + 1) * 
        teta)))))) * ((-1/teta) * (((1 - p1^-teta)^delta + (1 - 
    (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - p1^-teta)^(delta - 1) * (delta * (p1^-(teta + 1) * 
        teta)))))))

                  
 c.copula2.be2 <-(1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((1 - p1^-teta)^delta + (1 - 
    (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta)))))) * ((-1/teta) * (((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta)))))) + (1 - ((1 - p1^-teta)^delta + (1 - (1 - 
    p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 1) * ((-1/teta) * 
    (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - (1 - p2)^-teta)^(delta - 1) * 
        (delta * ((1 - p2)^-(teta + 1 + 1) * (teta + 1) * teta)) - 
        (1 - (1 - p2)^-teta)^((delta - 1) - 1) * ((delta - 1) * 
            ((1 - p2)^-(teta + 1) * teta)) * (delta * ((1 - p2)^-(teta + 
            1) * teta)))) - ((1 - p1^-teta)^delta + (1 - (1 - 
        p2)^-teta)^delta - 1)^(((1/delta) - 1) - 1) * (((1/delta) - 
        1) * ((1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - 
        p2)^-(teta + 1) * teta)))) * ((1/delta) * ((1 - (1 - 
        p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta))))))

                 



c.copula2.be1be2 <- -((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((1 - p1^-teta)^delta + (1 - 
    (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta)))))) * ((-1/teta) * (((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - p1^-teta)^(delta - 1) * (delta * (p1^-(teta + 1) * 
        teta)))))) - (1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^(1/delta))^((-1/teta) - 1) * ((-1/teta) * (((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^(((1/delta) - 1) - 1) * (((1/delta) - 
    1) * ((1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
    1) * teta)))) * ((1/delta) * ((1 - p1^-teta)^(delta - 1) * 
    (delta * (p1^-(teta + 1) * teta)))))))




c.copula2.be1th <-(-(((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^(1/delta))^((-1/teta) - 1) * (log((1 - ((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))) * (1/teta^2)) - 
    (1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
        1)^(1/delta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 1) * 
        (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
            1)^((1/delta) - 1) * ((1/delta) * ((1 - p1^-teta)^(delta - 
            1) * (delta * (p1^-teta * log(p1))) + (1 - (1 - p2)^-teta)^(delta - 
            1) * (delta * ((1 - p2)^-teta * log((1 - p2))))))))) * 
    ((-1/teta) * (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
        1)^((1/delta) - 1) * ((1/delta) * ((1 - p1^-teta)^(delta - 
        1) * (delta * (p1^-(teta + 1) * teta)))))) + (1 - ((1 - 
    p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 
    1) * (1/teta^2 * (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - p1^-teta)^(delta - 
    1) * (delta * (p1^-(teta + 1) * teta))))) + (-1/teta) * (((1 - 
    p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((1 - p1^-teta)^(delta - 1) * 
    (delta * (p1^-teta * log(p1))) + (1 - (1 - p2)^-teta)^(delta - 
    1) * (delta * ((1 - p2)^-teta * log((1 - p2)))))) * ((1/delta) * 
    ((1 - p1^-teta)^(delta - 1) * (delta * (p1^-(teta + 1) * 
        teta)))) + ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - p1^-teta)^((delta - 
    1) - 1) * ((delta - 1) * (p1^-teta * log(p1))) * (delta * 
    (p1^-(teta + 1) * teta)) + (1 - p1^-teta)^(delta - 1) * (delta * 
    (p1^-(teta + 1) - p1^-(teta + 1) * log(p1) * teta)))))))
)*(-exp(teta.st))


        


c.copula2.be2th <-(((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 
    1) * (log((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^(1/delta))) * (1/teta^2)) - (1 - ((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((1 - p1^-teta)^delta + (1 - 
    (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - p1^-teta)^(delta - 1) * (delta * (p1^-teta * log(p1))) + 
        (1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-teta * 
            log((1 - p2))))))))) * ((-1/teta) * (((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
        1) * teta)))))) + (1 - ((1 - p1^-teta)^delta + (1 - (1 - 
    p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 1) * (1/teta^2 * 
    (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - (1 - p2)^-teta)^(delta - 1) * 
        (delta * ((1 - p2)^-(teta + 1) * teta))))) + (-1/teta) * 
    (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((1 - p1^-teta)^(delta - 
        1) * (delta * (p1^-teta * log(p1))) + (1 - (1 - p2)^-teta)^(delta - 
        1) * (delta * ((1 - p2)^-teta * log((1 - p2)))))) * ((1/delta) * 
        ((1 - (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
            1) * teta)))) + ((1 - p1^-teta)^delta + (1 - (1 - 
        p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
        ((1 - (1 - p2)^-teta)^((delta - 1) - 1) * ((delta - 1) * 
            ((1 - p2)^-teta * log((1 - p2)))) * (delta * ((1 - 
            p2)^-(teta + 1) * teta)) + (1 - (1 - p2)^-teta)^(delta - 
            1) * (delta * ((1 - p2)^-(teta + 1) - (1 - p2)^-(teta + 
            1) * log((1 - p2)) * teta)))))))*(-exp(teta.st))





bit1.th2 <-((1 - ((1 - p1^(exp(teta.st) + 1))^delta + (1 - (1 - p2)^(exp(teta.st) + 
    1))^delta - 1)^(1/delta))^(((1/(exp(teta.st) + 1)) - 1) - 
    1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - p1^(exp(teta.st) + 
    1))^delta + (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^((1/delta) - 
    1) * ((1/delta) * ((1 - (1 - p2)^(exp(teta.st) + 1))^(delta - 
    1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
    exp(teta.st)))) + (1 - p1^(exp(teta.st) + 1))^(delta - 1) * 
    (delta * (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))))))) - 
    (1 - ((1 - p1^(exp(teta.st) + 1))^delta + (1 - (1 - p2)^(exp(teta.st) + 
        1))^delta - 1)^(1/delta))^((1/(exp(teta.st) + 1)) - 1) * 
        (log((1 - ((1 - p1^(exp(teta.st) + 1))^delta + (1 - (1 - 
            p2)^(exp(teta.st) + 1))^delta - 1)^(1/delta))) * 
            (exp(teta.st)/(exp(teta.st) + 1)^2))) * ((1/(exp(teta.st) + 
    1)) * (((1 - p1^(exp(teta.st) + 1))^delta + (1 - (1 - p2)^(exp(teta.st) + 
    1))^delta - 1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - 
    p2)^(exp(teta.st) + 1))^(delta - 1) * (delta * ((1 - p2)^(exp(teta.st) + 
    1) * (log((1 - p2)) * exp(teta.st)))) + (1 - p1^(exp(teta.st) + 
    1))^(delta - 1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * 
    exp(teta.st)))))))) + (1 - ((1 - p1^(exp(teta.st) + 1))^delta + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^(1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - p1^(exp(teta.st) + 
    1))^delta + (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^((1/delta) - 
    1) * ((1/delta) * ((1 - (1 - p2)^(exp(teta.st) + 1))^(delta - 
    1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
    exp(teta.st)) * (log((1 - p2)) * exp(teta.st)) + (1 - p2)^(exp(teta.st) + 
    1) * (log((1 - p2)) * exp(teta.st)))) - (1 - (1 - p2)^(exp(teta.st) + 
    1))^((delta - 1) - 1) * ((delta - 1) * ((1 - p2)^(exp(teta.st) + 
    1) * (log((1 - p2)) * exp(teta.st)))) * (delta * ((1 - p2)^(exp(teta.st) + 
    1) * (log((1 - p2)) * exp(teta.st)))) + ((1 - p1^(exp(teta.st) + 
    1))^(delta - 1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * 
    exp(teta.st)) * (log(p1) * exp(teta.st)) + p1^(exp(teta.st) + 
    1) * (log(p1) * exp(teta.st)))) - (1 - p1^(exp(teta.st) + 
    1))^((delta - 1) - 1) * ((delta - 1) * (p1^(exp(teta.st) + 
    1) * (log(p1) * exp(teta.st)))) * (delta * (p1^(exp(teta.st) + 
    1) * (log(p1) * exp(teta.st))))))) - ((1 - p1^(exp(teta.st) + 
    1))^delta + (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((1 - (1 - p2)^(exp(teta.st) + 
    1))^(delta - 1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * 
    (log((1 - p2)) * exp(teta.st)))) + (1 - p1^(exp(teta.st) + 
    1))^(delta - 1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * 
    exp(teta.st)))))) * ((1/delta) * ((1 - (1 - p2)^(exp(teta.st) + 
    1))^(delta - 1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * 
    (log((1 - p2)) * exp(teta.st)))) + (1 - p1^(exp(teta.st) + 
    1))^(delta - 1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * 
    exp(teta.st))))))) - exp(teta.st)/(exp(teta.st) + 1)^2 * 
    (((1 - p1^(exp(teta.st) + 1))^delta + (1 - (1 - p2)^(exp(teta.st) + 
        1))^delta - 1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - 
        p2)^(exp(teta.st) + 1))^(delta - 1) * (delta * ((1 - 
        p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)))) + 
        (1 - p1^(exp(teta.st) + 1))^(delta - 1) * (delta * (p1^(exp(teta.st) + 
            1) * (log(p1) * exp(teta.st)))))))) - (((1 - ((1 - 
    p1^(exp(teta.st) + 1))^delta + (1 - (1 - p2)^(exp(teta.st) + 
    1))^delta - 1)^(1/delta))^((1/(exp(teta.st) + 1)) - 1) * 
    ((1/(exp(teta.st) + 1)) * (((1 - p1^(exp(teta.st) + 1))^delta + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - (1 - p2)^(exp(teta.st) + 1))^(delta - 
        1) * (delta * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
        p2)) * exp(teta.st)))) + (1 - p1^(exp(teta.st) + 1))^(delta - 
        1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))))))) - 
    (1 - ((1 - p1^(exp(teta.st) + 1))^delta + (1 - (1 - p2)^(exp(teta.st) + 
        1))^delta - 1)^(1/delta))^(1/(exp(teta.st) + 1)) * (log((1 - 
        ((1 - p1^(exp(teta.st) + 1))^delta + (1 - (1 - p2)^(exp(teta.st) + 
            1))^delta - 1)^(1/delta))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2))) * (log((1 - ((1 - p1^(exp(teta.st) + 1))^delta + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^(1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) + (1 - ((1 - p1^(exp(teta.st) + 
    1))^delta + (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^(1/delta))^(1/(exp(teta.st) + 
    1)) * (((1 - p1^(exp(teta.st) + 1))^delta + (1 - (1 - p2)^(exp(teta.st) + 
    1))^delta - 1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - 
    p2)^(exp(teta.st) + 1))^(delta - 1) * (delta * ((1 - p2)^(exp(teta.st) + 
    1) * (log((1 - p2)) * exp(teta.st)))) + (1 - p1^(exp(teta.st) + 
    1))^(delta - 1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * 
    exp(teta.st))))))/(1 - ((1 - p1^(exp(teta.st) + 1))^delta + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^(1/delta)) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2) + log((1 - ((1 - p1^(exp(teta.st) + 
    1))^delta + (1 - (1 - p2)^(exp(teta.st) + 1))^delta - 1)^(1/delta))) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * (2 * 
        (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
        1)^2)^2)))





bit1.del2 <--((1 - ((1 - p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - 
    p2)^-teta)^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
    epsilon)))^((-1/teta) - 1) * ((-1/teta) * ((((1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) - 
    1)^(1/-(exp(delta.st) + epsilon)) * (log(((1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) - 
    1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - ((1 - 
    p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
    epsilon) - 1)^((1/-(exp(delta.st) + epsilon)) - 1) * ((1/-(exp(delta.st) + 
    epsilon)) * ((1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) * 
    (log((1 - (1 - p2)^-teta)) * exp(delta.st)) + (1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^-teta)) * exp(delta.st))))) * (log(((1 - 
    p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) + 
    ((1 - p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
        epsilon) - 1)^(1/-(exp(delta.st) + epsilon)) * (log(((1 - 
        p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
        epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2 - 
        exp(delta.st) * (2 * (exp(delta.st) * (exp(delta.st) + 
            epsilon)))/((-(exp(delta.st) + epsilon))^2)^2) - 
        ((1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) * (log((1 - 
            (1 - p2)^-teta)) * exp(delta.st)) + (1 - p1^-teta)^-(exp(delta.st) + 
            epsilon) * (log((1 - p1^-teta)) * exp(delta.st)))/((1 - 
            p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - 
            p2)^-teta)^-(exp(delta.st) + epsilon) - 1) * (exp(delta.st)/(-(exp(delta.st) + 
            epsilon))^2)) - ((((1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) - 
    1)^((1/-(exp(delta.st) + epsilon)) - 1) * (log(((1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) - 
    1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - ((1 - 
    p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
    epsilon) - 1)^(((1/-(exp(delta.st) + epsilon)) - 1) - 1) * 
    (((1/-(exp(delta.st) + epsilon)) - 1) * ((1 - (1 - p2)^-teta)^-(exp(delta.st) + 
        epsilon) * (log((1 - (1 - p2)^-teta)) * exp(delta.st)) + 
        (1 - p1^-teta)^-(exp(delta.st) + epsilon) * (log((1 - 
            p1^-teta)) * exp(delta.st))))) * ((1/-(exp(delta.st) + 
    epsilon)) * ((1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) * 
    (log((1 - (1 - p2)^-teta)) * exp(delta.st)) + (1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^-teta)) * exp(delta.st)))) + ((1 - 
    p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
    epsilon) - 1)^((1/-(exp(delta.st) + epsilon)) - 1) * (exp(delta.st)/(-(exp(delta.st) + 
    epsilon))^2 * ((1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) * 
    (log((1 - (1 - p2)^-teta)) * exp(delta.st)) + (1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^-teta)) * exp(delta.st))) + (1/-(exp(delta.st) + 
    epsilon)) * ((1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) * 
    (log((1 - (1 - p2)^-teta)) * exp(delta.st)) - (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - (1 - p2)^-teta)) * exp(delta.st)) * 
    (log((1 - (1 - p2)^-teta)) * exp(delta.st)) + ((1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^-teta)) * exp(delta.st)) - (1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^-teta)) * exp(delta.st)) * (log((1 - 
    p1^-teta)) * exp(delta.st)))))))) - (1 - ((1 - p1^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) - 
    1)^(1/-(exp(delta.st) + epsilon)))^(((-1/teta) - 1) - 1) * 
    (((-1/teta) - 1) * (((1 - p1^-teta)^-(exp(delta.st) + epsilon) + 
        (1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
        epsilon)) * (log(((1 - p1^-teta)^-(exp(delta.st) + epsilon) + 
        (1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) - 1)) * 
        (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - ((1 - 
        p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
        epsilon) - 1)^((1/-(exp(delta.st) + epsilon)) - 1) * 
        ((1/-(exp(delta.st) + epsilon)) * ((1 - (1 - p2)^-teta)^-(exp(delta.st) + 
            epsilon) * (log((1 - (1 - p2)^-teta)) * exp(delta.st)) + 
            (1 - p1^-teta)^-(exp(delta.st) + epsilon) * (log((1 - 
                p1^-teta)) * exp(delta.st)))))) * ((-1/teta) * 
    (((1 - p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
        epsilon) - 1)^(1/-(exp(delta.st) + epsilon)) * (log(((1 - 
        p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^-teta)^-(exp(delta.st) + 
        epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
        ((1 - p1^-teta)^-(exp(delta.st) + epsilon) + (1 - (1 - 
            p2)^-teta)^-(exp(delta.st) + epsilon) - 1)^((1/-(exp(delta.st) + 
            epsilon)) - 1) * ((1/-(exp(delta.st) + epsilon)) * 
            ((1 - (1 - p2)^-teta)^-(exp(delta.st) + epsilon) * 
                (log((1 - (1 - p2)^-teta)) * exp(delta.st)) + 
                (1 - p1^-teta)^-(exp(delta.st) + epsilon) * (log((1 - 
                  p1^-teta)) * exp(delta.st)))))))



c.copula2.be1del <-(-((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * ((((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((1 - p1^-teta)^delta * 
    log((1 - p1^-teta)) + (1 - (1 - p2)^-teta)^delta * log((1 - 
    (1 - p2)^-teta)))) - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^((1/delta) - 1) * (log(((1 - p1^-teta)^delta + (1 - (1 - 
    p2)^-teta)^delta - 1)) * (1/delta^2))) * ((1/delta) * ((1 - 
    p1^-teta)^(delta - 1) * (delta * (p1^-(teta + 1) * teta)))) + 
    ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - p1^-teta)^(delta - 1) * log((1 - 
        p1^-teta)) * (delta * (p1^-(teta + 1) * teta)) + (1 - 
        p1^-teta)^(delta - 1) * (p1^-(teta + 1) * teta)) - 1/delta^2 * 
        ((1 - p1^-teta)^(delta - 1) * (delta * (p1^-(teta + 1) * 
            teta)))))) - (1 - ((1 - p1^-teta)^delta + (1 - (1 - 
    p2)^-teta)^delta - 1)^(1/delta))^(((-1/teta) - 1) - 1) * 
    (((-1/teta) - 1) * (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
        1)^((1/delta) - 1) * ((1/delta) * ((1 - p1^-teta)^delta * 
        log((1 - p1^-teta)) + (1 - (1 - p2)^-teta)^delta * log((1 - 
        (1 - p2)^-teta)))) - ((1 - p1^-teta)^delta + (1 - (1 - 
        p2)^-teta)^delta - 1)^(1/delta) * (log(((1 - p1^-teta)^delta + 
        (1 - (1 - p2)^-teta)^delta - 1)) * (1/delta^2)))) * ((-1/teta) * 
    (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - p1^-teta)^(delta - 1) * (delta * 
        (p1^-(teta + 1) * teta))))))))*(-exp(delta.st))



        
c.copula2.be2del <- ((1 - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^((-1/teta) - 
    1) * ((-1/teta) * ((((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((1 - p1^-teta)^delta * 
    log((1 - p1^-teta)) + (1 - (1 - p2)^-teta)^delta * log((1 - 
    (1 - p2)^-teta)))) - ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^((1/delta) - 1) * (log(((1 - p1^-teta)^delta + (1 - (1 - 
    p2)^-teta)^delta - 1)) * (1/delta^2))) * ((1/delta) * ((1 - 
    (1 - p2)^-teta)^(delta - 1) * (delta * ((1 - p2)^-(teta + 
    1) * teta)))) + ((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - p2)^-teta)^(delta - 
    1) * log((1 - (1 - p2)^-teta)) * (delta * ((1 - p2)^-(teta + 
    1) * teta)) + (1 - (1 - p2)^-teta)^(delta - 1) * ((1 - p2)^-(teta + 
    1) * teta)) - 1/delta^2 * ((1 - (1 - p2)^-teta)^(delta - 
    1) * (delta * ((1 - p2)^-(teta + 1) * teta)))))) - (1 - ((1 - 
    p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^(1/delta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((1 - p1^-teta)^delta + (1 - 
    (1 - p2)^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - p1^-teta)^delta * log((1 - p1^-teta)) + (1 - (1 - p2)^-teta)^delta * 
        log((1 - (1 - p2)^-teta)))) - ((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)^(1/delta) * (log(((1 - p1^-teta)^delta + 
    (1 - (1 - p2)^-teta)^delta - 1)) * (1/delta^2)))) * ((-1/teta) * 
    (((1 - p1^-teta)^delta + (1 - (1 - p2)^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - (1 - p2)^-teta)^(delta - 1) * 
        (delta * ((1 - p2)^-(teta + 1) * teta)))))))*(-exp(delta.st))



bit1.thdel <--((1 - ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
    (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
    1)^(1/-(exp(delta.st) + epsilon)))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * ((((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^((1/-(exp(delta.st) + 
    epsilon)) - 1) * (log(((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^(((1/-(exp(delta.st) + epsilon)) - 1) - 1) * (((1/-(exp(delta.st) + 
        epsilon)) - 1) * ((1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - (1 - p2)^(exp(teta.st) + 1))) * 
        exp(delta.st)) + (1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st))))) * 
    ((1/-(exp(delta.st) + epsilon)) * ((1 - p1^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
        epsilon) * (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-((exp(delta.st) + 
            epsilon) + 1) * ((exp(delta.st) + epsilon) * ((1 - 
            p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)))))) + 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^((1/-(exp(delta.st) + epsilon)) - 1) * (exp(delta.st)/(-(exp(delta.st) + 
        epsilon))^2 * ((1 - p1^(exp(teta.st) + 1))^-((exp(delta.st) + 
        epsilon) + 1) * ((exp(delta.st) + epsilon) * (p1^(exp(teta.st) + 
        1) * (log(p1) * exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
        epsilon) * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * 
        exp(teta.st))))) + (1/-(exp(delta.st) + epsilon)) * ((1 - 
        p1^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 
        1) * (exp(delta.st) * (p1^(exp(teta.st) + 1) * (log(p1) * 
        exp(teta.st)))) - (1 - p1^(exp(teta.st) + 1))^-((exp(delta.st) + 
        epsilon) + 1) * (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st)) * 
        ((exp(delta.st) + epsilon) * (p1^(exp(teta.st) + 1) * 
            (log(p1) * exp(teta.st)))) + ((1 - (1 - p2)^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * (exp(delta.st) * 
        ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st)))) - 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-((exp(delta.st) + 
            epsilon) + 1) * (log((1 - (1 - p2)^(exp(teta.st) + 
            1))) * exp(delta.st)) * ((exp(delta.st) + epsilon) * 
            ((1 - p2)^(exp(teta.st) + 1) * (log((1 - p2)) * exp(teta.st))))))))) - 
    (1 - ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^(1/-(exp(delta.st) + epsilon)))^(((1/(exp(teta.st) + 
        1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - 
        p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + (1 - 
        (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^(1/-(exp(delta.st) + epsilon)) * (log(((1 - p1^(exp(teta.st) + 
        1))^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^(exp(teta.st) + 
        1))^-(exp(delta.st) + epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + 
        epsilon))^2)) - ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) - 1)^((1/-(exp(delta.st) + epsilon)) - 1) * 
        ((1/-(exp(delta.st) + epsilon)) * ((1 - (1 - p2)^(exp(teta.st) + 
            1))^-(exp(delta.st) + epsilon) * (log((1 - (1 - p2)^(exp(teta.st) + 
            1))) * exp(delta.st)) + (1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) * (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st)))))) * 
        ((1/(exp(teta.st) + 1)) * (((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) - 1)^((1/-(exp(delta.st) + epsilon)) - 1) * 
            ((1/-(exp(delta.st) + epsilon)) * ((1 - p1^(exp(teta.st) + 
                1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
                epsilon) * (p1^(exp(teta.st) + 1) * (log(p1) * 
                exp(teta.st)))) + (1 - (1 - p2)^(exp(teta.st) + 
                1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
                epsilon) * ((1 - p2)^(exp(teta.st) + 1) * (log((1 - 
                p2)) * exp(teta.st)))))))) - ((1 - ((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
    epsilon)))^(1/(exp(teta.st) + 1)) * ((((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
    epsilon)) * (log(((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^((1/-(exp(delta.st) + epsilon)) - 1) * ((1/-(exp(delta.st) + 
        epsilon)) * ((1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - (1 - p2)^(exp(teta.st) + 1))) * 
        exp(delta.st)) + (1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st)))))/(1 - 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^(1/-(exp(delta.st) + epsilon))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)^(1/-(exp(delta.st) + epsilon)))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - (1 - p2)^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
    epsilon)) * (log(((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^((1/-(exp(delta.st) + epsilon)) - 1) * ((1/-(exp(delta.st) + 
        epsilon)) * ((1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - (1 - p2)^(exp(teta.st) + 1))) * 
        exp(delta.st)) + (1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st)))))) * 
    (log((1 - ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) + (1 - (1 - p2)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) - 1)^(1/-(exp(delta.st) + epsilon)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2))))


}


if(BivD=="BB7.180"){
   
  c.copula.be1 <- 1 + (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^((-1/delta) - 1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 
    1) * (delta * (p1^(teta - 1) * teta))))))


 
  c.copula.be2 <- 1 + (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^((-1/delta) - 1) * ((-1/delta) * ((1 - p2^teta)^-(delta + 
    1) * (delta * (p2^(teta - 1) * teta))))))




  c.copula.theta <-((1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^(1/teta) * 
    (log((1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
        1)^(-1/delta))) * (1/teta^2)) + (1 - ((1 - p1^teta)^-delta + 
    (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 1) * ((1/teta) * 
    (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 1) * (delta * 
        (p1^teta * log(p1))) + (1 - p2^teta)^-(delta + 1) * (delta * 
        (p2^teta * log(p2))))))))*exp(teta.st)


  c.copula.delta <- ((1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(-1/delta) * (log(((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)) * (1/delta^2)) - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^((-1/delta) - 1) * ((-1/delta) * ((1 - p2^teta)^-delta * 
    log((1 - p2^teta)) + (1 - p1^teta)^-delta * log((1 - p1^teta)))))))*exp(delta.st)



 c.copula2.be1 <- (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((1 - p1^teta)^-(delta + 
    1) * (delta * (p1^(teta - 1) * teta)))) * ((-1/delta) * ((1 - 
    p1^teta)^-(delta + 1) * (delta * (p1^(teta - 1) * teta)))) + 
    ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 1 + 1) * 
        ((delta + 1) * (p1^(teta - 1) * teta)) * (delta * (p1^(teta - 
        1) * teta)) + (1 - p1^teta)^-(delta + 1) * (delta * (p1^((teta - 
        1) - 1) * (teta - 1) * teta)))))) - (1 - ((1 - p1^teta)^-delta + 
    (1 - p2^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 1) - 1) * 
    (((1/teta) - 1) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 
        1) * (delta * (p1^(teta - 1) * teta)))))) * ((1/teta) * 
    (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 1) * (delta * 
        (p1^(teta - 1) * teta))))))


                  
 c.copula2.be2 <- (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((1 - p2^teta)^-(delta + 
    1) * (delta * (p2^(teta - 1) * teta)))) * ((-1/delta) * ((1 - 
    p2^teta)^-(delta + 1) * (delta * (p2^(teta - 1) * teta)))) + 
    ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p2^teta)^-(delta + 1 + 1) * 
        ((delta + 1) * (p2^(teta - 1) * teta)) * (delta * (p2^(teta - 
        1) * teta)) + (1 - p2^teta)^-(delta + 1) * (delta * (p2^((teta - 
        1) - 1) * (teta - 1) * teta)))))) - (1 - ((1 - p1^teta)^-delta + 
    (1 - p2^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 1) - 1) * 
    (((1/teta) - 1) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
        1)^((-1/delta) - 1) * ((-1/delta) * ((1 - p2^teta)^-(delta + 
        1) * (delta * (p2^(teta - 1) * teta)))))) * ((1/teta) * 
    (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p2^teta)^-(delta + 1) * (delta * 
        (p2^(teta - 1) * teta))))))
                  
                  
                 




c.copula2.be1be2 <- (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((1 - p2^teta)^-(delta + 
    1) * (delta * (p2^(teta - 1) * teta)))) * ((-1/delta) * ((1 - 
    p1^teta)^-(delta + 1) * (delta * (p1^(teta - 1) * teta)))))) - 
    (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * (((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
        ((1 - p2^teta)^-(delta + 1) * (delta * (p2^(teta - 1) * 
            teta)))))) * ((1/teta) * (((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
        ((1 - p1^teta)^-(delta + 1) * (delta * (p1^(teta - 1) * 
            teta))))))





c.copula2.be1th <-((1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((1 - p1^teta)^-(delta + 
    1) * (delta * (p1^teta * log(p1))) + (1 - p2^teta)^-(delta + 
    1) * (delta * (p2^teta * log(p2))))) * ((-1/delta) * ((1 - 
    p1^teta)^-(delta + 1) * (delta * (p1^(teta - 1) * teta)))) + 
    ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 1 + 1) * 
        ((delta + 1) * (p1^teta * log(p1))) * (delta * (p1^(teta - 
        1) * teta)) + (1 - p1^teta)^-(delta + 1) * (delta * (p1^(teta - 
        1) * log(p1) * teta + p1^(teta - 1)))))) - 1/teta^2 * 
    (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 1) * (delta * 
        (p1^(teta - 1) * teta)))))) - ((1 - ((1 - p1^teta)^-delta + 
    (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 1) * (log((1 - 
    ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))) * 
    (1/teta^2)) + (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 1) * (delta * 
    (p1^teta * log(p1))) + (1 - p2^teta)^-(delta + 1) * (delta * 
    (p2^teta * log(p2)))))))) * ((1/teta) * (((1 - p1^teta)^-delta + 
    (1 - p2^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
    ((1 - p1^teta)^-(delta + 1) * (delta * (p1^(teta - 1) * teta)))))))*exp(teta.st)


        


c.copula2.be2th <-((1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((1 - p1^teta)^-(delta + 
    1) * (delta * (p1^teta * log(p1))) + (1 - p2^teta)^-(delta + 
    1) * (delta * (p2^teta * log(p2))))) * ((-1/delta) * ((1 - 
    p2^teta)^-(delta + 1) * (delta * (p2^(teta - 1) * teta)))) + 
    ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p2^teta)^-(delta + 1 + 1) * 
        ((delta + 1) * (p2^teta * log(p2))) * (delta * (p2^(teta - 
        1) * teta)) + (1 - p2^teta)^-(delta + 1) * (delta * (p2^(teta - 
        1) * log(p2) * teta + p2^(teta - 1)))))) - 1/teta^2 * 
    (((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p2^teta)^-(delta + 1) * (delta * 
        (p2^(teta - 1) * teta)))))) - ((1 - ((1 - p1^teta)^-delta + 
    (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 1) * (log((1 - 
    ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))) * 
    (1/teta^2)) + (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(-1/delta))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - p1^teta)^-(delta + 1) * (delta * 
    (p1^teta * log(p1))) + (1 - p2^teta)^-(delta + 1) * (delta * 
    (p2^teta * log(p2)))))))) * ((1/teta) * (((1 - p1^teta)^-delta + 
    (1 - p2^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
    ((1 - p2^teta)^-(delta + 1) * (delta * (p2^(teta - 1) * teta)))))))*exp(teta.st)






bit1.th2 <-(1 - ((1 - p1^(exp(teta.st) + 1))^-delta + (1 - p2^(exp(teta.st) + 
    1))^-delta - 1)^(-1/delta))^(1/(exp(teta.st) + 1)) * (log((1 - 
    ((1 - p1^(exp(teta.st) + 1))^-delta + (1 - p2^(exp(teta.st) + 
        1))^-delta - 1)^(-1/delta))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
    1)))/((exp(teta.st) + 1)^2)^2) - ((1 - p1^(exp(teta.st) + 
    1))^-delta + (1 - p2^(exp(teta.st) + 1))^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - p1^(exp(teta.st) + 1))^-(delta + 
    1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))) + 
    (1 - p2^(exp(teta.st) + 1))^-(delta + 1) * (delta * (p2^(exp(teta.st) + 
        1) * (log(p2) * exp(teta.st))))))/(1 - ((1 - p1^(exp(teta.st) + 
    1))^-delta + (1 - p2^(exp(teta.st) + 1))^-delta - 1)^(-1/delta)) * 
    (exp(teta.st)/(exp(teta.st) + 1)^2)) - ((1 - ((1 - p1^(exp(teta.st) + 
    1))^-delta + (1 - p2^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))^(1/(exp(teta.st) + 
    1)) * (log((1 - ((1 - p1^(exp(teta.st) + 1))^-delta + (1 - 
    p2^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - ((1 - p1^(exp(teta.st) + 1))^-delta + (1 - 
    p2^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - p1^(exp(teta.st) + 
    1))^-delta + (1 - p2^(exp(teta.st) + 1))^-delta - 1)^((-1/delta) - 
    1) * ((-1/delta) * ((1 - p1^(exp(teta.st) + 1))^-(delta + 
    1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))) + 
    (1 - p2^(exp(teta.st) + 1))^-(delta + 1) * (delta * (p2^(exp(teta.st) + 
        1) * (log(p2) * exp(teta.st))))))))) * (log((1 - ((1 - 
    p1^(exp(teta.st) + 1))^-delta + (1 - p2^(exp(teta.st) + 1))^-delta - 
    1)^(-1/delta))) * (exp(teta.st)/(exp(teta.st) + 1)^2)) + 
    ((1 - ((1 - p1^(exp(teta.st) + 1))^-delta + (1 - p2^(exp(teta.st) + 
        1))^-delta - 1)^(-1/delta))^((1/(exp(teta.st) + 1)) - 
        1) * ((1/(exp(teta.st) + 1)) * (((1 - p1^(exp(teta.st) + 
        1))^-delta + (1 - p2^(exp(teta.st) + 1))^-delta - 1)^(((-1/delta) - 
        1) - 1) * (((-1/delta) - 1) * ((1 - p1^(exp(teta.st) + 
        1))^-(delta + 1) * (delta * (p1^(exp(teta.st) + 1) * 
        (log(p1) * exp(teta.st)))) + (1 - p2^(exp(teta.st) + 
        1))^-(delta + 1) * (delta * (p2^(exp(teta.st) + 1) * 
        (log(p2) * exp(teta.st)))))) * ((-1/delta) * ((1 - p1^(exp(teta.st) + 
        1))^-(delta + 1) * (delta * (p1^(exp(teta.st) + 1) * 
        (log(p1) * exp(teta.st)))) + (1 - p2^(exp(teta.st) + 
        1))^-(delta + 1) * (delta * (p2^(exp(teta.st) + 1) * 
        (log(p2) * exp(teta.st)))))) + ((1 - p1^(exp(teta.st) + 
        1))^-delta + (1 - p2^(exp(teta.st) + 1))^-delta - 1)^((-1/delta) - 
        1) * ((-1/delta) * ((1 - p1^(exp(teta.st) + 1))^-(delta + 
        1 + 1) * ((delta + 1) * (p1^(exp(teta.st) + 1) * (log(p1) * 
        exp(teta.st)))) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * 
        exp(teta.st)))) + (1 - p1^(exp(teta.st) + 1))^-(delta + 
        1) * (delta * (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)) * 
        (log(p1) * exp(teta.st)) + p1^(exp(teta.st) + 1) * (log(p1) * 
        exp(teta.st)))) + ((1 - p2^(exp(teta.st) + 1))^-(delta + 
        1 + 1) * ((delta + 1) * (p2^(exp(teta.st) + 1) * (log(p2) * 
        exp(teta.st)))) * (delta * (p2^(exp(teta.st) + 1) * (log(p2) * 
        exp(teta.st)))) + (1 - p2^(exp(teta.st) + 1))^-(delta + 
        1) * (delta * (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)) * 
        (log(p2) * exp(teta.st)) + p2^(exp(teta.st) + 1) * (log(p2) * 
        exp(teta.st)))))))) - exp(teta.st)/(exp(teta.st) + 1)^2 * 
        (((1 - p1^(exp(teta.st) + 1))^-delta + (1 - p2^(exp(teta.st) + 
            1))^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
            ((1 - p1^(exp(teta.st) + 1))^-(delta + 1) * (delta * 
                (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))) + 
                (1 - p2^(exp(teta.st) + 1))^-(delta + 1) * (delta * 
                  (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)))))))) - 
        ((1 - ((1 - p1^(exp(teta.st) + 1))^-delta + (1 - p2^(exp(teta.st) + 
            1))^-delta - 1)^(-1/delta))^((1/(exp(teta.st) + 1)) - 
            1) * (log((1 - ((1 - p1^(exp(teta.st) + 1))^-delta + 
            (1 - p2^(exp(teta.st) + 1))^-delta - 1)^(-1/delta))) * 
            (exp(teta.st)/(exp(teta.st) + 1)^2)) + (1 - ((1 - 
            p1^(exp(teta.st) + 1))^-delta + (1 - p2^(exp(teta.st) + 
            1))^-delta - 1)^(-1/delta))^(((1/(exp(teta.st) + 
            1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - 
            p1^(exp(teta.st) + 1))^-delta + (1 - p2^(exp(teta.st) + 
            1))^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
            ((1 - p1^(exp(teta.st) + 1))^-(delta + 1) * (delta * 
                (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))) + 
                (1 - p2^(exp(teta.st) + 1))^-(delta + 1) * (delta * 
                  (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st))))))))) * 
            ((1/(exp(teta.st) + 1)) * (((1 - p1^(exp(teta.st) + 
                1))^-delta + (1 - p2^(exp(teta.st) + 1))^-delta - 
                1)^((-1/delta) - 1) * ((-1/delta) * ((1 - p1^(exp(teta.st) + 
                1))^-(delta + 1) * (delta * (p1^(exp(teta.st) + 
                1) * (log(p1) * exp(teta.st)))) + (1 - p2^(exp(teta.st) + 
                1))^-(delta + 1) * (delta * (p2^(exp(teta.st) + 
                1) * (log(p2) * exp(teta.st)))))))))




bit1.del2 <-(1 - ((1 - p1^teta)^-(exp(delta.st) + epsilon) + (1 - p2^teta)^-(exp(delta.st) + 
    epsilon) - 1)^(-1/(exp(delta.st) + epsilon)))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1^teta)^-(exp(delta.st) + epsilon) + 
    (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)^(-1/(exp(delta.st) + 
    epsilon)) * (log(((1 - p1^teta)^-(exp(delta.st) + epsilon) + 
    (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)) * (exp(delta.st)/(exp(delta.st) + 
    epsilon)^2)) - ((1 - p1^teta)^-(exp(delta.st) + epsilon) + 
    (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)^((-1/(exp(delta.st) + 
    epsilon)) - 1) * ((-1/(exp(delta.st) + epsilon)) * ((1 - 
    p2^teta)^-(exp(delta.st) + epsilon) * (log((1 - p2^teta)) * 
    exp(delta.st)) + (1 - p1^teta)^-(exp(delta.st) + epsilon) * 
    (log((1 - p1^teta)) * exp(delta.st))))) * (log(((1 - p1^teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)) * 
    (exp(delta.st)/(exp(delta.st) + epsilon)^2)) + ((1 - p1^teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)^(-1/(exp(delta.st) + 
    epsilon)) * (log(((1 - p1^teta)^-(exp(delta.st) + epsilon) + 
    (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)) * (exp(delta.st)/(exp(delta.st) + 
    epsilon)^2 - exp(delta.st) * (2 * (exp(delta.st) * (exp(delta.st) + 
    epsilon)))/((exp(delta.st) + epsilon)^2)^2) - ((1 - p2^teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p2^teta)) * exp(delta.st)) + (1 - p1^teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^teta)) * exp(delta.st)))/((1 - p1^teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1) * 
    (exp(delta.st)/(exp(delta.st) + epsilon)^2)) - ((((1 - p1^teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)^((-1/(exp(delta.st) + 
    epsilon)) - 1) * (log(((1 - p1^teta)^-(exp(delta.st) + epsilon) + 
    (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)) * (exp(delta.st)/(exp(delta.st) + 
    epsilon)^2)) - ((1 - p1^teta)^-(exp(delta.st) + epsilon) + 
    (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)^(((-1/(exp(delta.st) + 
    epsilon)) - 1) - 1) * (((-1/(exp(delta.st) + epsilon)) - 
    1) * ((1 - p2^teta)^-(exp(delta.st) + epsilon) * (log((1 - 
    p2^teta)) * exp(delta.st)) + (1 - p1^teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^teta)) * exp(delta.st))))) * ((-1/(exp(delta.st) + 
    epsilon)) * ((1 - p2^teta)^-(exp(delta.st) + epsilon) * (log((1 - 
    p2^teta)) * exp(delta.st)) + (1 - p1^teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^teta)) * exp(delta.st)))) + ((1 - 
    p1^teta)^-(exp(delta.st) + epsilon) + (1 - p2^teta)^-(exp(delta.st) + 
    epsilon) - 1)^((-1/(exp(delta.st) + epsilon)) - 1) * (exp(delta.st)/(exp(delta.st) + 
    epsilon)^2 * ((1 - p2^teta)^-(exp(delta.st) + epsilon) * 
    (log((1 - p2^teta)) * exp(delta.st)) + (1 - p1^teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^teta)) * exp(delta.st))) + (-1/(exp(delta.st) + 
    epsilon)) * ((1 - p2^teta)^-(exp(delta.st) + epsilon) * (log((1 - 
    p2^teta)) * exp(delta.st)) - (1 - p2^teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p2^teta)) * exp(delta.st)) * (log((1 - 
    p2^teta)) * exp(delta.st)) + ((1 - p1^teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^teta)) * exp(delta.st)) - (1 - p1^teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - p1^teta)) * exp(delta.st)) * (log((1 - 
    p1^teta)) * exp(delta.st)))))))) - (1 - ((1 - p1^teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^teta)^-(exp(delta.st) + epsilon) - 1)^(-1/(exp(delta.st) + 
    epsilon)))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    p1^teta)^-(exp(delta.st) + epsilon) + (1 - p2^teta)^-(exp(delta.st) + 
    epsilon) - 1)^(-1/(exp(delta.st) + epsilon)) * (log(((1 - 
    p1^teta)^-(exp(delta.st) + epsilon) + (1 - p2^teta)^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(exp(delta.st) + epsilon)^2)) - 
    ((1 - p1^teta)^-(exp(delta.st) + epsilon) + (1 - p2^teta)^-(exp(delta.st) + 
        epsilon) - 1)^((-1/(exp(delta.st) + epsilon)) - 1) * 
        ((-1/(exp(delta.st) + epsilon)) * ((1 - p2^teta)^-(exp(delta.st) + 
            epsilon) * (log((1 - p2^teta)) * exp(delta.st)) + 
            (1 - p1^teta)^-(exp(delta.st) + epsilon) * (log((1 - 
                p1^teta)) * exp(delta.st)))))) * ((1/teta) * 
    (((1 - p1^teta)^-(exp(delta.st) + epsilon) + (1 - p2^teta)^-(exp(delta.st) + 
        epsilon) - 1)^(-1/(exp(delta.st) + epsilon)) * (log(((1 - 
        p1^teta)^-(exp(delta.st) + epsilon) + (1 - p2^teta)^-(exp(delta.st) + 
        epsilon) - 1)) * (exp(delta.st)/(exp(delta.st) + epsilon)^2)) - 
        ((1 - p1^teta)^-(exp(delta.st) + epsilon) + (1 - p2^teta)^-(exp(delta.st) + 
            epsilon) - 1)^((-1/(exp(delta.st) + epsilon)) - 1) * 
            ((-1/(exp(delta.st) + epsilon)) * ((1 - p2^teta)^-(exp(delta.st) + 
                epsilon) * (log((1 - p2^teta)) * exp(delta.st)) + 
                (1 - p1^teta)^-(exp(delta.st) + epsilon) * (log((1 - 
                  p1^teta)) * exp(delta.st))))))



c.copula2.be1del <- ((1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^((-1/delta) - 1) * (log(((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)) * (1/delta^2)) - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((1 - p2^teta)^-delta * 
    log((1 - p2^teta)) + (1 - p1^teta)^-delta * log((1 - p1^teta))))) * 
    ((-1/delta) * ((1 - p1^teta)^-(delta + 1) * (delta * (p1^(teta - 
        1) * teta)))) + ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^((-1/delta) - 1) * (1/delta^2 * ((1 - p1^teta)^-(delta + 
    1) * (delta * (p1^(teta - 1) * teta))) + (-1/delta) * ((1 - 
    p1^teta)^-(delta + 1) * (p1^(teta - 1) * teta) - (1 - p1^teta)^-(delta + 
    1) * log((1 - p1^teta)) * (delta * (p1^(teta - 1) * teta)))))) - 
    (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * (((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)^(-1/delta) * (log(((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)) * (1/delta^2)) - ((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
        ((1 - p2^teta)^-delta * log((1 - p2^teta)) + (1 - p1^teta)^-delta * 
            log((1 - p1^teta)))))) * ((1/teta) * (((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
        ((1 - p1^teta)^-(delta + 1) * (delta * (p1^(teta - 1) * 
            teta)))))))*exp(delta.st)


        



        
c.copula2.be2del <- ((1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^((1/teta) - 
    1) * ((1/teta) * ((((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^((-1/delta) - 1) * (log(((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)) * (1/delta^2)) - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^(((-1/delta) - 1) - 1) * (((-1/delta) - 1) * ((1 - p2^teta)^-delta * 
    log((1 - p2^teta)) + (1 - p1^teta)^-delta * log((1 - p1^teta))))) * 
    ((-1/delta) * ((1 - p2^teta)^-(delta + 1) * (delta * (p2^(teta - 
        1) * teta)))) + ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 
    1)^((-1/delta) - 1) * (1/delta^2 * ((1 - p2^teta)^-(delta + 
    1) * (delta * (p2^(teta - 1) * teta))) + (-1/delta) * ((1 - 
    p2^teta)^-(delta + 1) * (p2^(teta - 1) * teta) - (1 - p2^teta)^-(delta + 
    1) * log((1 - p2^teta)) * (delta * (p2^(teta - 1) * teta)))))) - 
    (1 - ((1 - p1^teta)^-delta + (1 - p2^teta)^-delta - 1)^(-1/delta))^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * (((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)^(-1/delta) * (log(((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)) * (1/delta^2)) - ((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
        ((1 - p2^teta)^-delta * log((1 - p2^teta)) + (1 - p1^teta)^-delta * 
            log((1 - p1^teta)))))) * ((1/teta) * (((1 - p1^teta)^-delta + 
        (1 - p2^teta)^-delta - 1)^((-1/delta) - 1) * ((-1/delta) * 
        ((1 - p2^teta)^-(delta + 1) * (delta * (p2^(teta - 1) * 
            teta)))))))*exp(delta.st)



bit1.thdel <-(1 - ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
    (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
    1)^(-1/(exp(delta.st) + epsilon)))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * ((((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^((-1/(exp(delta.st) + 
    epsilon)) - 1) * (log(((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(exp(delta.st) + epsilon)^2)) - 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^(((-1/(exp(delta.st) + epsilon)) - 1) - 1) * (((-1/(exp(delta.st) + 
        epsilon)) - 1) * ((1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - p2^(exp(teta.st) + 1))) * exp(delta.st)) + 
        (1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) * 
            (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st))))) * 
    ((-1/(exp(delta.st) + epsilon)) * ((1 - p1^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
        epsilon) * (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))) + 
        (1 - p2^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 
            1) * ((exp(delta.st) + epsilon) * (p2^(exp(teta.st) + 
            1) * (log(p2) * exp(teta.st)))))) + ((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^((-1/(exp(delta.st) + 
    epsilon)) - 1) * (exp(delta.st)/(exp(delta.st) + epsilon)^2 * 
    ((1 - p1^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 
        1) * ((exp(delta.st) + epsilon) * (p1^(exp(teta.st) + 
        1) * (log(p1) * exp(teta.st)))) + (1 - p2^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
        epsilon) * (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st))))) + 
    (-1/(exp(delta.st) + epsilon)) * ((1 - p1^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * (exp(delta.st) * 
        (p1^(exp(teta.st) + 1) * (log(p1) * exp(teta.st)))) - 
        (1 - p1^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 
            1) * (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st)) * 
            ((exp(delta.st) + epsilon) * (p1^(exp(teta.st) + 
                1) * (log(p1) * exp(teta.st)))) + ((1 - p2^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * (exp(delta.st) * 
        (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)))) - 
        (1 - p2^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 
            1) * (log((1 - p2^(exp(teta.st) + 1))) * exp(delta.st)) * 
            ((exp(delta.st) + epsilon) * (p2^(exp(teta.st) + 
                1) * (log(p2) * exp(teta.st))))))))) - (1 - ((1 - 
    p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + (1 - 
    p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 1)^(-1/(exp(delta.st) + 
    epsilon)))^(((1/(exp(teta.st) + 1)) - 1) - 1) * (((1/(exp(teta.st) + 
    1)) - 1) * (((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)^(-1/(exp(delta.st) + epsilon)) * (log(((1 - 
    p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + (1 - 
    p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 1)) * 
    (exp(delta.st)/(exp(delta.st) + epsilon)^2)) - ((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^((-1/(exp(delta.st) + 
    epsilon)) - 1) * ((-1/(exp(delta.st) + epsilon)) * ((1 - 
    p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) * (log((1 - 
    p2^(exp(teta.st) + 1))) * exp(delta.st)) + (1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) * (log((1 - p1^(exp(teta.st) + 
    1))) * exp(delta.st)))))) * ((1/(exp(teta.st) + 1)) * (((1 - 
    p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + (1 - 
    p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 1)^((-1/(exp(delta.st) + 
    epsilon)) - 1) * ((-1/(exp(delta.st) + epsilon)) * ((1 - 
    p1^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 1) * 
    ((exp(delta.st) + epsilon) * (p1^(exp(teta.st) + 1) * (log(p1) * 
        exp(teta.st)))) + (1 - p2^(exp(teta.st) + 1))^-((exp(delta.st) + 
    epsilon) + 1) * ((exp(delta.st) + epsilon) * (p2^(exp(teta.st) + 
    1) * (log(p2) * exp(teta.st)))))))) - ((1 - ((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^(-1/(exp(delta.st) + 
    epsilon)))^(1/(exp(teta.st) + 1)) * ((((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^(-1/(exp(delta.st) + 
    epsilon)) * (log(((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(exp(delta.st) + epsilon)^2)) - 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^((-1/(exp(delta.st) + epsilon)) - 1) * ((-1/(exp(delta.st) + 
        epsilon)) * ((1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - p2^(exp(teta.st) + 1))) * exp(delta.st)) + 
        (1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) * 
            (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st)))))/(1 - 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^(-1/(exp(delta.st) + epsilon))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)^(-1/(exp(delta.st) + epsilon)))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - p1^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^(-1/(exp(delta.st) + 
    epsilon)) * (log(((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(exp(delta.st) + epsilon)^2)) - 
    ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^((-1/(exp(delta.st) + epsilon)) - 1) * ((-1/(exp(delta.st) + 
        epsilon)) * ((1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - p2^(exp(teta.st) + 1))) * exp(delta.st)) + 
        (1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) * 
            (log((1 - p1^(exp(teta.st) + 1))) * exp(delta.st)))))) * 
    (log((1 - ((1 - p1^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) - 1)^(-1/(exp(delta.st) + epsilon)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)))

}


if(BivD=="BB7.270"){
  
   
  c.copula.be1 <- 1 + (1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(1/delta))^((1/-teta) - 1) * ((1/-teta) * (((1 - (1 - 
    p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^((1/delta) - 
    1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 1) * (delta * 
    ((1 - p1)^-(teta + 1) * teta))))))



 
  c.copula.be2 <- -((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 
    1) * ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - p2^-teta)^(delta - 
    1) * (delta * (p2^-(teta + 1) * teta)))))))


  c.copula.theta <- ((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^(1/-teta) * 
    (log((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
        1)^(1/delta))) * (1/(-teta)^2)) - (1 - ((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 1) * ((1/-teta) * 
    (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 1) * 
        (delta * ((1 - p1)^-teta * log((1 - p1)))) + (1 - p2^-teta)^(delta - 
        1) * (delta * (p2^-teta * log(p2))))))))*(-exp(teta.st))


  c.copula.delta <- (-((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 
    1) * ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - p1)^-teta)^delta * 
    log((1 - (1 - p1)^-teta)) + (1 - p2^-teta)^delta * log((1 - 
    p2^-teta)))) - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(1/delta) * (log(((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)) * (1/delta^2))))))*(-exp(delta.st))


 c.copula2.be1 <-  (1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^(((1/-teta) - 
    1) - 1) * (((1/-teta) - 1) * (((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p1)^-teta)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta)))))) * ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p1)^-teta)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta)))))) + (1 - ((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 1) * ((1/-teta) * 
    (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 1) * 
        (delta * ((1 - p1)^-(teta + 1 + 1) * (teta + 1) * teta)) - 
        (1 - (1 - p1)^-teta)^((delta - 1) - 1) * ((delta - 1) * 
            ((1 - p1)^-(teta + 1) * teta)) * (delta * ((1 - p1)^-(teta + 
            1) * teta)))) - ((1 - (1 - p1)^-teta)^delta + (1 - 
        p2^-teta)^delta - 1)^(((1/delta) - 1) - 1) * (((1/delta) - 
        1) * ((1 - (1 - p1)^-teta)^(delta - 1) * (delta * ((1 - 
        p1)^-(teta + 1) * teta)))) * ((1/delta) * ((1 - (1 - 
        p1)^-teta)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta))))))


                  
 c.copula2.be2 <- -((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 
    1) * ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((1 - p2^-teta)^(delta - 
    1) * (delta * (p2^-(teta + 1) * teta)))) * ((1/delta) * ((1 - 
    p2^-teta)^(delta - 1) * (delta * (p2^-(teta + 1) * teta)))) + 
    ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - p2^-teta)^((delta - 1) - 1) * 
        ((delta - 1) * (p2^-(teta + 1) * teta)) * (delta * (p2^-(teta + 
        1) * teta)) - (1 - p2^-teta)^(delta - 1) * (delta * (p2^-(teta + 
        1 + 1) * (teta + 1) * teta)))))) - (1 - ((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^(1/delta))^(((1/-teta) - 1) - 1) * 
    (((1/-teta) - 1) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
        1)^((1/delta) - 1) * ((1/delta) * ((1 - p2^-teta)^(delta - 
        1) * (delta * (p2^-(teta + 1) * teta)))))) * ((1/-teta) * 
    (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - p2^-teta)^(delta - 1) * (delta * 
        (p2^-(teta + 1) * teta)))))))

                  
                 



c.copula2.be1be2 <- (1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 
    1) * ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((1 - p2^-teta)^(delta - 
    1) * (delta * (p2^-(teta + 1) * teta)))) * ((1/delta) * ((1 - 
    (1 - p1)^-teta)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
    1) * teta)))))) - (1 - ((1 - (1 - p1)^-teta)^delta + (1 - 
    p2^-teta)^delta - 1)^(1/delta))^(((1/-teta) - 1) - 1) * (((1/-teta) - 
    1) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - p2^-teta)^(delta - 
    1) * (delta * (p2^-(teta + 1) * teta)))))) * ((1/-teta) * 
    (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 1) * 
        (delta * ((1 - p1)^-(teta + 1) * teta))))))





c.copula2.be1th <-(((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 
    1) * (log((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(1/delta))) * (1/(-teta)^2)) - (1 - ((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^(1/delta))^(((1/-teta) - 1) - 1) * 
    (((1/-teta) - 1) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
        1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 
        1) * (delta * ((1 - p1)^-teta * log((1 - p1)))) + (1 - 
        p2^-teta)^(delta - 1) * (delta * (p2^-teta * log(p2)))))))) * 
    ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
        1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 
        1) * (delta * ((1 - p1)^-(teta + 1) * teta)))))) + (1 - 
    ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 
    1) * (1/(-teta)^2 * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 
    1) * (delta * ((1 - p1)^-(teta + 1) * teta))))) + (1/-teta) * 
    (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(((1/delta) - 
        1) - 1) * (((1/delta) - 1) * ((1 - (1 - p1)^-teta)^(delta - 
        1) * (delta * ((1 - p1)^-teta * log((1 - p1)))) + (1 - 
        p2^-teta)^(delta - 1) * (delta * (p2^-teta * log(p2))))) * 
        ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 1) * (delta * 
            ((1 - p1)^-(teta + 1) * teta)))) + ((1 - (1 - p1)^-teta)^delta + 
        (1 - p2^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
        ((1 - (1 - p1)^-teta)^((delta - 1) - 1) * ((delta - 1) * 
            ((1 - p1)^-teta * log((1 - p1)))) * (delta * ((1 - 
            p1)^-(teta + 1) * teta)) + (1 - (1 - p1)^-teta)^(delta - 
            1) * (delta * ((1 - p1)^-(teta + 1) - (1 - p1)^-(teta + 
            1) * log((1 - p1)) * teta)))))))*(-exp(teta.st))


        


c.copula2.be2th <-(-(((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(1/delta))^((1/-teta) - 1) * (log((1 - ((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^(1/delta))) * (1/(-teta)^2)) - 
    (1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
        1)^(1/delta))^(((1/-teta) - 1) - 1) * (((1/-teta) - 1) * 
        (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
            1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 
            1) * (delta * ((1 - p1)^-teta * log((1 - p1)))) + 
            (1 - p2^-teta)^(delta - 1) * (delta * (p2^-teta * 
                log(p2)))))))) * ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - p2^-teta)^(delta - 1) * (delta * (p2^-(teta + 1) * 
        teta)))))) + (1 - ((1 - (1 - p1)^-teta)^delta + (1 - 
    p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 1) * (1/(-teta)^2 * 
    (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - p2^-teta)^(delta - 1) * (delta * 
        (p2^-(teta + 1) * teta))))) + (1/-teta) * (((1 - (1 - 
    p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((1 - (1 - p1)^-teta)^(delta - 
    1) * (delta * ((1 - p1)^-teta * log((1 - p1)))) + (1 - p2^-teta)^(delta - 
    1) * (delta * (p2^-teta * log(p2))))) * ((1/delta) * ((1 - 
    p2^-teta)^(delta - 1) * (delta * (p2^-(teta + 1) * teta)))) + 
    ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - p2^-teta)^((delta - 1) - 1) * 
        ((delta - 1) * (p2^-teta * log(p2))) * (delta * (p2^-(teta + 
        1) * teta)) + (1 - p2^-teta)^(delta - 1) * (delta * (p2^-(teta + 
        1) - p2^-(teta + 1) * log(p2) * teta))))))))*(-exp(teta.st))





bit1.th2 <-((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^delta + (1 - p2^(exp(teta.st) + 
    1))^delta - 1)^(1/delta))^(((1/(exp(teta.st) + 1)) - 1) - 
    1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 - p1)^(exp(teta.st) + 
    1))^delta + (1 - p2^(exp(teta.st) + 1))^delta - 1)^((1/delta) - 
    1) * ((1/delta) * ((1 - p2^(exp(teta.st) + 1))^(delta - 1) * 
    (delta * (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)))) + 
    (1 - (1 - p1)^(exp(teta.st) + 1))^(delta - 1) * (delta * 
        ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st)))))))) - 
    (1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^delta + (1 - p2^(exp(teta.st) + 
        1))^delta - 1)^(1/delta))^((1/(exp(teta.st) + 1)) - 1) * 
        (log((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^delta + 
            (1 - p2^(exp(teta.st) + 1))^delta - 1)^(1/delta))) * 
            (exp(teta.st)/(exp(teta.st) + 1)^2))) * ((1/(exp(teta.st) + 
    1)) * (((1 - (1 - p1)^(exp(teta.st) + 1))^delta + (1 - p2^(exp(teta.st) + 
    1))^delta - 1)^((1/delta) - 1) * ((1/delta) * ((1 - p2^(exp(teta.st) + 
    1))^(delta - 1) * (delta * (p2^(exp(teta.st) + 1) * (log(p2) * 
    exp(teta.st)))) + (1 - (1 - p1)^(exp(teta.st) + 1))^(delta - 
    1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
    exp(teta.st)))))))) + (1 - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^delta + (1 - p2^(exp(teta.st) + 1))^delta - 1)^(1/delta))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - p1)^(exp(teta.st) + 
    1))^delta + (1 - p2^(exp(teta.st) + 1))^delta - 1)^((1/delta) - 
    1) * ((1/delta) * ((1 - p2^(exp(teta.st) + 1))^(delta - 1) * 
    (delta * (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)) * 
        (log(p2) * exp(teta.st)) + p2^(exp(teta.st) + 1) * (log(p2) * 
        exp(teta.st)))) - (1 - p2^(exp(teta.st) + 1))^((delta - 
    1) - 1) * ((delta - 1) * (p2^(exp(teta.st) + 1) * (log(p2) * 
    exp(teta.st)))) * (delta * (p2^(exp(teta.st) + 1) * (log(p2) * 
    exp(teta.st)))) + ((1 - (1 - p1)^(exp(teta.st) + 1))^(delta - 
    1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
    exp(teta.st)) * (log((1 - p1)) * exp(teta.st)) + (1 - p1)^(exp(teta.st) + 
    1) * (log((1 - p1)) * exp(teta.st)))) - (1 - (1 - p1)^(exp(teta.st) + 
    1))^((delta - 1) - 1) * ((delta - 1) * ((1 - p1)^(exp(teta.st) + 
    1) * (log((1 - p1)) * exp(teta.st)))) * (delta * ((1 - p1)^(exp(teta.st) + 
    1) * (log((1 - p1)) * exp(teta.st))))))) - ((1 - (1 - p1)^(exp(teta.st) + 
    1))^delta + (1 - p2^(exp(teta.st) + 1))^delta - 1)^(((1/delta) - 
    1) - 1) * (((1/delta) - 1) * ((1 - p2^(exp(teta.st) + 1))^(delta - 
    1) * (delta * (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)))) + 
    (1 - (1 - p1)^(exp(teta.st) + 1))^(delta - 1) * (delta * 
        ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st)))))) * 
    ((1/delta) * ((1 - p2^(exp(teta.st) + 1))^(delta - 1) * (delta * 
        (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)))) + 
        (1 - (1 - p1)^(exp(teta.st) + 1))^(delta - 1) * (delta * 
            ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st))))))) - 
    exp(teta.st)/(exp(teta.st) + 1)^2 * (((1 - (1 - p1)^(exp(teta.st) + 
        1))^delta + (1 - p2^(exp(teta.st) + 1))^delta - 1)^((1/delta) - 
        1) * ((1/delta) * ((1 - p2^(exp(teta.st) + 1))^(delta - 
        1) * (delta * (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)))) + 
        (1 - (1 - p1)^(exp(teta.st) + 1))^(delta - 1) * (delta * 
            ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st)))))))) - 
    (((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^delta + (1 - p2^(exp(teta.st) + 
        1))^delta - 1)^(1/delta))^((1/(exp(teta.st) + 1)) - 1) * 
        ((1/(exp(teta.st) + 1)) * (((1 - (1 - p1)^(exp(teta.st) + 
            1))^delta + (1 - p2^(exp(teta.st) + 1))^delta - 1)^((1/delta) - 
            1) * ((1/delta) * ((1 - p2^(exp(teta.st) + 1))^(delta - 
            1) * (delta * (p2^(exp(teta.st) + 1) * (log(p2) * 
            exp(teta.st)))) + (1 - (1 - p1)^(exp(teta.st) + 1))^(delta - 
            1) * (delta * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - 
            p1)) * exp(teta.st)))))))) - (1 - ((1 - (1 - p1)^(exp(teta.st) + 
        1))^delta + (1 - p2^(exp(teta.st) + 1))^delta - 1)^(1/delta))^(1/(exp(teta.st) + 
        1)) * (log((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^delta + 
        (1 - p2^(exp(teta.st) + 1))^delta - 1)^(1/delta))) * 
        (exp(teta.st)/(exp(teta.st) + 1)^2))) * (log((1 - ((1 - 
        (1 - p1)^(exp(teta.st) + 1))^delta + (1 - p2^(exp(teta.st) + 
        1))^delta - 1)^(1/delta))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)) + (1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^delta + 
        (1 - p2^(exp(teta.st) + 1))^delta - 1)^(1/delta))^(1/(exp(teta.st) + 
        1)) * (((1 - (1 - p1)^(exp(teta.st) + 1))^delta + (1 - 
        p2^(exp(teta.st) + 1))^delta - 1)^((1/delta) - 1) * ((1/delta) * 
        ((1 - p2^(exp(teta.st) + 1))^(delta - 1) * (delta * (p2^(exp(teta.st) + 
            1) * (log(p2) * exp(teta.st)))) + (1 - (1 - p1)^(exp(teta.st) + 
            1))^(delta - 1) * (delta * ((1 - p1)^(exp(teta.st) + 
            1) * (log((1 - p1)) * exp(teta.st))))))/(1 - ((1 - 
        (1 - p1)^(exp(teta.st) + 1))^delta + (1 - p2^(exp(teta.st) + 
        1))^delta - 1)^(1/delta)) * (exp(teta.st)/(exp(teta.st) + 
        1)^2) + log((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^delta + 
        (1 - p2^(exp(teta.st) + 1))^delta - 1)^(1/delta))) * 
        (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * (2 * 
            (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
            1)^2)^2)))




bit1.del2 <--((1 - ((1 - (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - 
    p2^-teta)^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
    epsilon)))^((1/-teta) - 1) * ((1/-teta) * ((((1 - (1 - p1)^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
    epsilon)) * (log(((1 - (1 - p1)^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 1)) * 
    (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - ((1 - (1 - 
    p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - p2^-teta)^-(exp(delta.st) + 
    epsilon) - 1)^((1/-(exp(delta.st) + epsilon)) - 1) * ((1/-(exp(delta.st) + 
    epsilon)) * ((1 - p2^-teta)^-(exp(delta.st) + epsilon) * 
    (log((1 - p2^-teta)) * exp(delta.st)) + (1 - (1 - p1)^-teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - (1 - p1)^-teta)) * exp(delta.st))))) * 
    (log(((1 - (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - 
        p2^-teta)^-(exp(delta.st) + epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + 
        epsilon))^2)) + ((1 - (1 - p1)^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
    epsilon)) * (log(((1 - (1 - p1)^-teta)^-(exp(delta.st) + 
    epsilon) + (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 1)) * 
    (exp(delta.st)/(-(exp(delta.st) + epsilon))^2 - exp(delta.st) * 
        (2 * (exp(delta.st) * (exp(delta.st) + epsilon)))/((-(exp(delta.st) + 
        epsilon))^2)^2) - ((1 - p2^-teta)^-(exp(delta.st) + epsilon) * 
    (log((1 - p2^-teta)) * exp(delta.st)) + (1 - (1 - p1)^-teta)^-(exp(delta.st) + 
    epsilon) * (log((1 - (1 - p1)^-teta)) * exp(delta.st)))/((1 - 
    (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - p2^-teta)^-(exp(delta.st) + 
    epsilon) - 1) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
    ((((1 - (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - 
        p2^-teta)^-(exp(delta.st) + epsilon) - 1)^((1/-(exp(delta.st) + 
        epsilon)) - 1) * (log(((1 - (1 - p1)^-teta)^-(exp(delta.st) + 
        epsilon) + (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 
        1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
        ((1 - (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - 
            p2^-teta)^-(exp(delta.st) + epsilon) - 1)^(((1/-(exp(delta.st) + 
            epsilon)) - 1) - 1) * (((1/-(exp(delta.st) + epsilon)) - 
            1) * ((1 - p2^-teta)^-(exp(delta.st) + epsilon) * 
            (log((1 - p2^-teta)) * exp(delta.st)) + (1 - (1 - 
            p1)^-teta)^-(exp(delta.st) + epsilon) * (log((1 - 
            (1 - p1)^-teta)) * exp(delta.st))))) * ((1/-(exp(delta.st) + 
        epsilon)) * ((1 - p2^-teta)^-(exp(delta.st) + epsilon) * 
        (log((1 - p2^-teta)) * exp(delta.st)) + (1 - (1 - p1)^-teta)^-(exp(delta.st) + 
        epsilon) * (log((1 - (1 - p1)^-teta)) * exp(delta.st)))) + 
        ((1 - (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - 
            p2^-teta)^-(exp(delta.st) + epsilon) - 1)^((1/-(exp(delta.st) + 
            epsilon)) - 1) * (exp(delta.st)/(-(exp(delta.st) + 
            epsilon))^2 * ((1 - p2^-teta)^-(exp(delta.st) + epsilon) * 
            (log((1 - p2^-teta)) * exp(delta.st)) + (1 - (1 - 
            p1)^-teta)^-(exp(delta.st) + epsilon) * (log((1 - 
            (1 - p1)^-teta)) * exp(delta.st))) + (1/-(exp(delta.st) + 
            epsilon)) * ((1 - p2^-teta)^-(exp(delta.st) + epsilon) * 
            (log((1 - p2^-teta)) * exp(delta.st)) - (1 - p2^-teta)^-(exp(delta.st) + 
            epsilon) * (log((1 - p2^-teta)) * exp(delta.st)) * 
            (log((1 - p2^-teta)) * exp(delta.st)) + ((1 - (1 - 
            p1)^-teta)^-(exp(delta.st) + epsilon) * (log((1 - 
            (1 - p1)^-teta)) * exp(delta.st)) - (1 - (1 - p1)^-teta)^-(exp(delta.st) + 
            epsilon) * (log((1 - (1 - p1)^-teta)) * exp(delta.st)) * 
            (log((1 - (1 - p1)^-teta)) * exp(delta.st)))))))) - 
    (1 - ((1 - (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - 
        p2^-teta)^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
        epsilon)))^(((1/-teta) - 1) - 1) * (((1/-teta) - 1) * 
        (((1 - (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - 
            p2^-teta)^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
            epsilon)) * (log(((1 - (1 - p1)^-teta)^-(exp(delta.st) + 
            epsilon) + (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 
            1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
            ((1 - (1 - p1)^-teta)^-(exp(delta.st) + epsilon) + 
                (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 1)^((1/-(exp(delta.st) + 
                epsilon)) - 1) * ((1/-(exp(delta.st) + epsilon)) * 
                ((1 - p2^-teta)^-(exp(delta.st) + epsilon) * 
                  (log((1 - p2^-teta)) * exp(delta.st)) + (1 - 
                  (1 - p1)^-teta)^-(exp(delta.st) + epsilon) * 
                  (log((1 - (1 - p1)^-teta)) * exp(delta.st)))))) * 
        ((1/-teta) * (((1 - (1 - p1)^-teta)^-(exp(delta.st) + 
            epsilon) + (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 
            1)^(1/-(exp(delta.st) + epsilon)) * (log(((1 - (1 - 
            p1)^-teta)^-(exp(delta.st) + epsilon) + (1 - p2^-teta)^-(exp(delta.st) + 
            epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + 
            epsilon))^2)) - ((1 - (1 - p1)^-teta)^-(exp(delta.st) + 
            epsilon) + (1 - p2^-teta)^-(exp(delta.st) + epsilon) - 
            1)^((1/-(exp(delta.st) + epsilon)) - 1) * ((1/-(exp(delta.st) + 
            epsilon)) * ((1 - p2^-teta)^-(exp(delta.st) + epsilon) * 
            (log((1 - p2^-teta)) * exp(delta.st)) + (1 - (1 - 
            p1)^-teta)^-(exp(delta.st) + epsilon) * (log((1 - 
            (1 - p1)^-teta)) * exp(delta.st)))))))





c.copula2.be1del <-((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 
    1) * ((1/-teta) * ((((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((1 - (1 - 
    p1)^-teta)^delta * log((1 - (1 - p1)^-teta)) + (1 - p2^-teta)^delta * 
    log((1 - p2^-teta)))) - ((1 - (1 - p1)^-teta)^delta + (1 - 
    p2^-teta)^delta - 1)^((1/delta) - 1) * (log(((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)) * (1/delta^2))) * ((1/delta) * 
    ((1 - (1 - p1)^-teta)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta)))) + ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - p1)^-teta)^(delta - 
    1) * log((1 - (1 - p1)^-teta)) * (delta * ((1 - p1)^-(teta + 
    1) * teta)) + (1 - (1 - p1)^-teta)^(delta - 1) * ((1 - p1)^-(teta + 
    1) * teta)) - 1/delta^2 * ((1 - (1 - p1)^-teta)^(delta - 
    1) * (delta * ((1 - p1)^-(teta + 1) * teta)))))) - (1 - ((1 - 
    (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^(((1/-teta) - 
    1) - 1) * (((1/-teta) - 1) * (((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p1)^-teta)^delta * log((1 - (1 - p1)^-teta)) + 
        (1 - p2^-teta)^delta * log((1 - p2^-teta)))) - ((1 - 
    (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta) * 
    (log(((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
        1)) * (1/delta^2)))) * ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - (1 - p1)^-teta)^(delta - 1) * (delta * ((1 - p1)^-(teta + 
        1) * teta)))))))*(-exp(delta.st))


        



        
c.copula2.be2del <- (-((1 - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 1)^(1/delta))^((1/-teta) - 
    1) * ((1/-teta) * ((((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(((1/delta) - 1) - 1) * (((1/delta) - 1) * ((1 - (1 - 
    p1)^-teta)^delta * log((1 - (1 - p1)^-teta)) + (1 - p2^-teta)^delta * 
    log((1 - p2^-teta)))) - ((1 - (1 - p1)^-teta)^delta + (1 - 
    p2^-teta)^delta - 1)^((1/delta) - 1) * (log(((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)) * (1/delta^2))) * ((1/delta) * 
    ((1 - p2^-teta)^(delta - 1) * (delta * (p2^-(teta + 1) * 
        teta)))) + ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - p2^-teta)^(delta - 
    1) * log((1 - p2^-teta)) * (delta * (p2^-(teta + 1) * teta)) + 
    (1 - p2^-teta)^(delta - 1) * (p2^-(teta + 1) * teta)) - 1/delta^2 * 
    ((1 - p2^-teta)^(delta - 1) * (delta * (p2^-(teta + 1) * 
        teta)))))) - (1 - ((1 - (1 - p1)^-teta)^delta + (1 - 
    p2^-teta)^delta - 1)^(1/delta))^(((1/-teta) - 1) - 1) * (((1/-teta) - 
    1) * (((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^((1/delta) - 1) * ((1/delta) * ((1 - (1 - p1)^-teta)^delta * 
    log((1 - (1 - p1)^-teta)) + (1 - p2^-teta)^delta * log((1 - 
    p2^-teta)))) - ((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)^(1/delta) * (log(((1 - (1 - p1)^-teta)^delta + (1 - p2^-teta)^delta - 
    1)) * (1/delta^2)))) * ((1/-teta) * (((1 - (1 - p1)^-teta)^delta + 
    (1 - p2^-teta)^delta - 1)^((1/delta) - 1) * ((1/delta) * 
    ((1 - p2^-teta)^(delta - 1) * (delta * (p2^-(teta + 1) * 
        teta))))))))*(-exp(delta.st))




bit1.thdel <--((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
    (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
    1)^(1/-(exp(delta.st) + epsilon)))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * ((((1 - (1 - p1)^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^((1/-(exp(delta.st) + 
    epsilon)) - 1) * (log(((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^(((1/-(exp(delta.st) + epsilon)) - 1) - 1) * (((1/-(exp(delta.st) + 
        epsilon)) - 1) * ((1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - p2^(exp(teta.st) + 1))) * exp(delta.st)) + 
        (1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) * 
            (log((1 - (1 - p1)^(exp(teta.st) + 1))) * exp(delta.st))))) * 
    ((1/-(exp(delta.st) + epsilon)) * ((1 - (1 - p1)^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
        epsilon) * ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * 
        exp(teta.st)))) + (1 - p2^(exp(teta.st) + 1))^-((exp(delta.st) + 
        epsilon) + 1) * ((exp(delta.st) + epsilon) * (p2^(exp(teta.st) + 
        1) * (log(p2) * exp(teta.st)))))) + ((1 - (1 - p1)^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^((1/-(exp(delta.st) + 
    epsilon)) - 1) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2 * 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 
        1) * ((exp(delta.st) + epsilon) * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st)))) + (1 - p2^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
        epsilon) * (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st))))) + 
    (1/-(exp(delta.st) + epsilon)) * ((1 - (1 - p1)^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * (exp(delta.st) * 
        ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st)))) - 
        (1 - (1 - p1)^(exp(teta.st) + 1))^-((exp(delta.st) + 
            epsilon) + 1) * (log((1 - (1 - p1)^(exp(teta.st) + 
            1))) * exp(delta.st)) * ((exp(delta.st) + epsilon) * 
            ((1 - p1)^(exp(teta.st) + 1) * (log((1 - p1)) * exp(teta.st)))) + 
        ((1 - p2^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 
            1) * (exp(delta.st) * (p2^(exp(teta.st) + 1) * (log(p2) * 
            exp(teta.st)))) - (1 - p2^(exp(teta.st) + 1))^-((exp(delta.st) + 
            epsilon) + 1) * (log((1 - p2^(exp(teta.st) + 1))) * 
            exp(delta.st)) * ((exp(delta.st) + epsilon) * (p2^(exp(teta.st) + 
            1) * (log(p2) * exp(teta.st))))))))) - (1 - ((1 - 
    (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
    (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
    1)^(1/-(exp(delta.st) + epsilon)))^(((1/(exp(teta.st) + 1)) - 
    1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 - p1)^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
    1))^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
    epsilon)) * (log(((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
    epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
    ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)^((1/-(exp(delta.st) + epsilon)) - 1) * ((1/-(exp(delta.st) + 
        epsilon)) * ((1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) * (log((1 - p2^(exp(teta.st) + 1))) * exp(delta.st)) + 
        (1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) * 
            (log((1 - (1 - p1)^(exp(teta.st) + 1))) * exp(delta.st)))))) * 
    ((1/(exp(teta.st) + 1)) * (((1 - (1 - p1)^(exp(teta.st) + 
        1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
        1))^-(exp(delta.st) + epsilon) - 1)^((1/-(exp(delta.st) + 
        epsilon)) - 1) * ((1/-(exp(delta.st) + epsilon)) * ((1 - 
        (1 - p1)^(exp(teta.st) + 1))^-((exp(delta.st) + epsilon) + 
        1) * ((exp(delta.st) + epsilon) * ((1 - p1)^(exp(teta.st) + 
        1) * (log((1 - p1)) * exp(teta.st)))) + (1 - p2^(exp(teta.st) + 
        1))^-((exp(delta.st) + epsilon) + 1) * ((exp(delta.st) + 
        epsilon) * (p2^(exp(teta.st) + 1) * (log(p2) * exp(teta.st)))))))) - 
    ((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) - 1)^(1/-(exp(delta.st) + epsilon)))^(1/(exp(teta.st) + 
        1)) * ((((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) - 1)^(1/-(exp(delta.st) + epsilon)) * (log(((1 - 
        (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) + 
        (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + epsilon) - 
        1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
        ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) - 1)^((1/-(exp(delta.st) + epsilon)) - 1) * 
            ((1/-(exp(delta.st) + epsilon)) * ((1 - p2^(exp(teta.st) + 
                1))^-(exp(delta.st) + epsilon) * (log((1 - p2^(exp(teta.st) + 
                1))) * exp(delta.st)) + (1 - (1 - p1)^(exp(teta.st) + 
                1))^-(exp(delta.st) + epsilon) * (log((1 - (1 - 
                p1)^(exp(teta.st) + 1))) * exp(delta.st)))))/(1 - 
        ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) - 1)^(1/-(exp(delta.st) + epsilon))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)) + (1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) - 1)^(1/-(exp(delta.st) + epsilon)))^((1/(exp(teta.st) + 
        1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - p1)^(exp(teta.st) + 
        1))^-(exp(delta.st) + epsilon) + (1 - p2^(exp(teta.st) + 
        1))^-(exp(delta.st) + epsilon) - 1)^(1/-(exp(delta.st) + 
        epsilon)) * (log(((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
        epsilon) - 1)) * (exp(delta.st)/(-(exp(delta.st) + epsilon))^2)) - 
        ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) - 1)^((1/-(exp(delta.st) + epsilon)) - 1) * 
            ((1/-(exp(delta.st) + epsilon)) * ((1 - p2^(exp(teta.st) + 
                1))^-(exp(delta.st) + epsilon) * (log((1 - p2^(exp(teta.st) + 
                1))) * exp(delta.st)) + (1 - (1 - p1)^(exp(teta.st) + 
                1))^-(exp(delta.st) + epsilon) * (log((1 - (1 - 
                p1)^(exp(teta.st) + 1))) * exp(delta.st)))))) * 
        (log((1 - ((1 - (1 - p1)^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) + (1 - p2^(exp(teta.st) + 1))^-(exp(delta.st) + 
            epsilon) - 1)^(1/-(exp(delta.st) + epsilon)))) * 
            (exp(teta.st)/(exp(teta.st) + 1)^2))))

}



if(BivD=="BB8.0"){

epsilon <- 0
  
   
  c.copula.be1 <- 1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * ((1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 
        1) * (teta * delta)) * (1 - (1 - delta * p2)^teta))))

 
  c.copula.be2 <-1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * ((1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
        p1)^teta) * ((1 - delta * p2)^(teta - 1) * (teta * delta)))))




  c.copula.theta <-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^(1/teta) * (log((1 - 
    (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * 
        (1 - (1 - delta * p2)^teta))) * (1/teta^2)) + (1 - (1 - 
    (1 - delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * (1 - 
    (1 - delta * p2)^teta))^((1/teta) - 1) * ((1/teta) * (((1 - 
    (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^teta * log((1 - 
    delta))) * (1 - (1 - delta * p1)^teta) - (1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * p1)^teta * log((1 - delta * p1)))) * (1 - (1 - 
    delta * p2)^teta) - (1 - (1 - delta)^teta)^-1 * (1 - (1 - 
    delta * p1)^teta) * ((1 - delta * p2)^teta * log((1 - delta * 
    p2)))))))*exp(teta.st)



  c.copula.delta <- (1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * (((1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 
        1) * (teta * p1)) - (1 - (1 - delta)^teta)^-(1 + 1) * 
        ((1 - delta)^(teta - 1) * teta) * (1 - (1 - delta * p1)^teta)) * 
        (1 - (1 - delta * p2)^teta) + (1 - (1 - delta)^teta)^-1 * 
        (1 - (1 - delta * p1)^teta) * ((1 - delta * p2)^(teta - 
        1) * (teta * p2))))) - 1/delta^2 * (1 - (1 - (1 - (1 - 
    delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * (1 - (1 - 
    delta * p2)^teta))^(1/teta)))*dnorm(delta.st)




 c.copula2.be1 <- -(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * ((1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^((teta - 
        1) - 1) * ((teta - 1) * delta) * (teta * delta)) * (1 - 
        (1 - delta * p2)^teta))) + (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * (1 - (1 - delta * p2)^teta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * p1)^(teta - 1) * (teta * delta)) * (1 - (1 - 
    delta * p2)^teta))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * p1)^(teta - 1) * (teta * delta)) * (1 - (1 - 
    delta * p2)^teta)))))


                  
 c.copula2.be2 <- -(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * ((1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
        p1)^teta) * ((1 - delta * p2)^((teta - 1) - 1) * ((teta - 
        1) * delta) * (teta * delta)))) + (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * (1 - (1 - delta * p2)^teta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * ((1 - delta * p2)^(teta - 1) * 
    (teta * delta)))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * ((1 - delta * p2)^(teta - 1) * 
    (teta * delta))))))



                  
                  




c.copula2.be1be2 <-1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * ((1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 
        1) * (teta * delta)) * ((1 - delta * p2)^(teta - 1) * 
        (teta * delta)))) - (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * (1 - (1 - delta * p2)^teta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * ((1 - delta * p2)^(teta - 1) * 
    (teta * delta)))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * p1)^(teta - 1) * (teta * delta)) * (1 - (1 - 
    delta * p2)^teta))))





c.copula2.be1th <-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * (((1 - (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^teta * 
        log((1 - delta))) * ((1 - delta * p1)^(teta - 1) * (teta * 
        delta)) + (1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 
        1) * log((1 - delta * p1)) * (teta * delta) + (1 - delta * 
        p1)^(teta - 1) * delta)) * (1 - (1 - delta * p2)^teta) - 
        (1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 
            1) * (teta * delta)) * ((1 - delta * p2)^teta * log((1 - 
            delta * p2)))) - 1/teta^2 * ((1 - (1 - delta)^teta)^-1 * 
        ((1 - delta * p1)^(teta - 1) * (teta * delta)) * (1 - 
        (1 - delta * p2)^teta))) - ((1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 
    1) * (log((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))) * (1/teta^2)) + 
    (1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * 
        (1 - (1 - delta * p2)^teta))^(((1/teta) - 1) - 1) * (((1/teta) - 
        1) * (((1 - (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^teta * 
        log((1 - delta))) * (1 - (1 - delta * p1)^teta) - (1 - 
        (1 - delta)^teta)^-1 * ((1 - delta * p1)^teta * log((1 - 
        delta * p1)))) * (1 - (1 - delta * p2)^teta) - (1 - (1 - 
        delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * ((1 - 
        delta * p2)^teta * log((1 - delta * p2)))))) * ((1/teta) * 
    ((1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 1) * 
        (teta * delta)) * (1 - (1 - delta * p2)^teta)))))*exp(teta.st)



c.copula2.be2th <-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * (((1 - (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^teta * 
        log((1 - delta))) * (1 - (1 - delta * p1)^teta) - (1 - 
        (1 - delta)^teta)^-1 * ((1 - delta * p1)^teta * log((1 - 
        delta * p1)))) * ((1 - delta * p2)^(teta - 1) * (teta * 
        delta)) + (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
        p1)^teta) * ((1 - delta * p2)^(teta - 1) * log((1 - delta * 
        p2)) * (teta * delta) + (1 - delta * p2)^(teta - 1) * 
        delta)) - 1/teta^2 * ((1 - (1 - delta)^teta)^-1 * (1 - 
        (1 - delta * p1)^teta) * ((1 - delta * p2)^(teta - 1) * 
        (teta * delta)))) - ((1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 
    1) * (log((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))) * (1/teta^2)) + 
    (1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * 
        (1 - (1 - delta * p2)^teta))^(((1/teta) - 1) - 1) * (((1/teta) - 
        1) * (((1 - (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^teta * 
        log((1 - delta))) * (1 - (1 - delta * p1)^teta) - (1 - 
        (1 - delta)^teta)^-1 * ((1 - delta * p1)^teta * log((1 - 
        delta * p1)))) * (1 - (1 - delta * p2)^teta) - (1 - (1 - 
        delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * ((1 - 
        delta * p2)^teta * log((1 - delta * p2)))))) * ((1/teta) * 
    ((1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * 
        ((1 - delta * p2)^(teta - 1) * (teta * delta))))))*exp(teta.st)









bit1.th2 <-1/delta * ((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - 
    (1 - delta * p1)^(exp(teta.st) + 1)) * (1 - (1 - delta * 
    p2)^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 1)) * (log((1 - 
    (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - delta * 
        p1)^(exp(teta.st) + 1)) * (1 - (1 - delta * p2)^(exp(teta.st) + 
        1)))) * (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * 
    (2 * (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
    1)^2)^2) - (((1 - (1 - delta)^(exp(teta.st) + 1))^-(1 + 1) * 
    ((1 - delta)^(exp(teta.st) + 1) * (log((1 - delta)) * exp(teta.st))) * 
    (1 - (1 - delta * p1)^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * ((1 - delta * p1)^(exp(teta.st) + 1) * (log((1 - 
    delta * p1)) * exp(teta.st)))) * (1 - (1 - delta * p2)^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - 
    delta * p1)^(exp(teta.st) + 1)) * ((1 - delta * p2)^(exp(teta.st) + 
    1) * (log((1 - delta * p2)) * exp(teta.st))))/(1 - (1 - (1 - 
    delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - delta * p1)^(exp(teta.st) + 
    1)) * (1 - (1 - delta * p2)^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) - ((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - delta * p1)^(exp(teta.st) + 1)) * (1 - (1 - delta * 
    p2)^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 1)) * (log((1 - 
    (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - delta * 
        p1)^(exp(teta.st) + 1)) * (1 - (1 - delta * p2)^(exp(teta.st) + 
        1)))) * (exp(teta.st)/(exp(teta.st) + 1)^2)) + (1 - (1 - 
    (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - delta * p1)^(exp(teta.st) + 
    1)) * (1 - (1 - delta * p2)^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) * (1 - (1 - delta * p1)^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * 
    p1)^(exp(teta.st) + 1) * (log((1 - delta * p1)) * exp(teta.st)))) * 
    (1 - (1 - delta * p2)^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 - delta * p1)^(exp(teta.st) + 1)) * ((1 - 
    delta * p2)^(exp(teta.st) + 1) * (log((1 - delta * p2)) * 
    exp(teta.st)))))) * (log((1 - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 - delta * p1)^(exp(teta.st) + 1)) * (1 - 
    (1 - delta * p2)^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + ((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - delta * p1)^(exp(teta.st) + 1)) * (1 - (1 - delta * 
    p2)^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 
    1)) * ((((1 - (1 - delta)^(exp(teta.st) + 1))^-(1 + 1 + 1) * 
    ((1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - delta)) * 
        exp(teta.st)))) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) + (1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st)) * (log((1 - delta)) * exp(teta.st)) + 
    (1 - delta)^(exp(teta.st) + 1) * (log((1 - delta)) * exp(teta.st)))) * 
    (1 - (1 - delta * p1)^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) * ((1 - delta * p1)^(exp(teta.st) + 
    1) * (log((1 - delta * p1)) * exp(teta.st))) - ((1 - (1 - 
    delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 
    1) * (log((1 - delta)) * exp(teta.st))) * ((1 - delta * p1)^(exp(teta.st) + 
    1) * (log((1 - delta * p1)) * exp(teta.st))) + (1 - (1 - 
    delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * p1)^(exp(teta.st) + 
    1) * (log((1 - delta * p1)) * exp(teta.st)) * (log((1 - delta * 
    p1)) * exp(teta.st)) + (1 - delta * p1)^(exp(teta.st) + 1) * 
    (log((1 - delta * p1)) * exp(teta.st))))) * (1 - (1 - delta * 
    p2)^(exp(teta.st) + 1)) - ((1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) * (1 - (1 - delta * p1)^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * 
    p1)^(exp(teta.st) + 1) * (log((1 - delta * p1)) * exp(teta.st)))) * 
    ((1 - delta * p2)^(exp(teta.st) + 1) * (log((1 - delta * 
        p2)) * exp(teta.st))) - (((1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) * (1 - (1 - delta * p1)^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * 
    p1)^(exp(teta.st) + 1) * (log((1 - delta * p1)) * exp(teta.st)))) * 
    ((1 - delta * p2)^(exp(teta.st) + 1) * (log((1 - delta * 
        p2)) * exp(teta.st))) + (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 - delta * p1)^(exp(teta.st) + 1)) * ((1 - 
    delta * p2)^(exp(teta.st) + 1) * (log((1 - delta * p2)) * 
    exp(teta.st)) * (log((1 - delta * p2)) * exp(teta.st)) + 
    (1 - delta * p2)^(exp(teta.st) + 1) * (log((1 - delta * p2)) * 
        exp(teta.st))))) - exp(teta.st)/(exp(teta.st) + 1)^2 * 
    (((1 - (1 - delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 
        1) * (log((1 - delta)) * exp(teta.st))) * (1 - (1 - delta * 
        p1)^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
        1))^-1 * ((1 - delta * p1)^(exp(teta.st) + 1) * (log((1 - 
        delta * p1)) * exp(teta.st)))) * (1 - (1 - delta * p2)^(exp(teta.st) + 
        1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - 
        (1 - delta * p1)^(exp(teta.st) + 1)) * ((1 - delta * 
        p2)^(exp(teta.st) + 1) * (log((1 - delta * p2)) * exp(teta.st))))) - 
    ((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - 
        delta * p1)^(exp(teta.st) + 1)) * (1 - (1 - delta * p2)^(exp(teta.st) + 
        1)))^((1/(exp(teta.st) + 1)) - 1) * (log((1 - (1 - (1 - 
        delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - delta * p1)^(exp(teta.st) + 
        1)) * (1 - (1 - delta * p2)^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)) + (1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
        (1 - (1 - delta * p1)^(exp(teta.st) + 1)) * (1 - (1 - 
        delta * p2)^(exp(teta.st) + 1)))^(((1/(exp(teta.st) + 
        1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - 
        (1 - delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 
        1) * (log((1 - delta)) * exp(teta.st))) * (1 - (1 - delta * 
        p1)^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
        1))^-1 * ((1 - delta * p1)^(exp(teta.st) + 1) * (log((1 - 
        delta * p1)) * exp(teta.st)))) * (1 - (1 - delta * p2)^(exp(teta.st) + 
        1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - 
        (1 - delta * p1)^(exp(teta.st) + 1)) * ((1 - delta * 
        p2)^(exp(teta.st) + 1) * (log((1 - delta * p2)) * exp(teta.st)))))) * 
        ((1/(exp(teta.st) + 1)) * (((1 - (1 - delta)^(exp(teta.st) + 
            1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * 
            (log((1 - delta)) * exp(teta.st))) * (1 - (1 - delta * 
            p1)^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
            1))^-1 * ((1 - delta * p1)^(exp(teta.st) + 1) * (log((1 - 
            delta * p1)) * exp(teta.st)))) * (1 - (1 - delta * 
            p2)^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
            1))^-1 * (1 - (1 - delta * p1)^(exp(teta.st) + 1)) * 
            ((1 - delta * p2)^(exp(teta.st) + 1) * (log((1 - 
                delta * p2)) * exp(teta.st)))))))



bit1.del2 <-1/(pnorm(delta.st) + epsilon) * ((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    p1)^teta) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^teta))^((1/teta) - 
    1) * ((1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
    ((1 - (pnorm(delta.st) + epsilon) * p1)^(teta - 1) * (teta * 
        (dnorm(delta.st) * p1))) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 
    1) * (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p1)^teta)) * ((1 - (pnorm(delta.st) + epsilon) * 
    p2)^(teta - 1) * (teta * (dnorm(delta.st) * p2))) - ((1 - 
    (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * ((1 - (pnorm(delta.st) + 
    epsilon) * p1)^(teta - 1) * (teta * (delta.st * dnorm(delta.st) * 
    p1)) + (1 - (pnorm(delta.st) + epsilon) * p1)^((teta - 1) - 
    1) * ((teta - 1) * (dnorm(delta.st) * p1)) * (teta * (dnorm(delta.st) * 
    p1))) + (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 
    1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * 
    dnorm(delta.st))) * ((1 - (pnorm(delta.st) + epsilon) * p1)^(teta - 
    1) * (teta * (dnorm(delta.st) * p1))) + ((1 - (1 - (pnorm(delta.st) + 
    epsilon))^teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 
    1) * (teta * dnorm(delta.st))) * ((1 - (pnorm(delta.st) + 
    epsilon) * p1)^(teta - 1) * (teta * (dnorm(delta.st) * p1))) - 
    ((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 1) * 
        ((1 - (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * 
            (delta.st * dnorm(delta.st))) + (1 - (pnorm(delta.st) + 
            epsilon))^((teta - 1) - 1) * ((teta - 1) * dnorm(delta.st)) * 
            (teta * dnorm(delta.st))) + (1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-(1 + 1 + 1) * ((1 + 1) * ((1 - (pnorm(delta.st) + 
        epsilon))^(teta - 1) * (teta * dnorm(delta.st)))) * ((1 - 
        (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * dnorm(delta.st)))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^teta))) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^teta) + (((1 - 
    (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * ((1 - (pnorm(delta.st) + 
    epsilon) * p1)^(teta - 1) * (teta * (dnorm(delta.st) * p1))) - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 1) * ((1 - 
        (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * dnorm(delta.st))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^teta)) * 
    ((1 - (pnorm(delta.st) + epsilon) * p2)^(teta - 1) * (teta * 
        (dnorm(delta.st) * p2))) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    p1)^teta) * ((1 - (pnorm(delta.st) + epsilon) * p2)^(teta - 
    1) * (teta * (delta.st * dnorm(delta.st) * p2)) + (1 - (pnorm(delta.st) + 
    epsilon) * p2)^((teta - 1) - 1) * ((teta - 1) * (dnorm(delta.st) * 
    p2)) * (teta * (dnorm(delta.st) * p2)))))) - (1 - (1 - (1 - 
    (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p1)^teta) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    p2)^teta))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * ((1 - (pnorm(delta.st) + 
    epsilon) * p1)^(teta - 1) * (teta * (dnorm(delta.st) * p1))) - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 1) * ((1 - 
        (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * dnorm(delta.st))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^teta)) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^teta) + (1 - 
    (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p1)^teta) * ((1 - (pnorm(delta.st) + epsilon) * 
    p2)^(teta - 1) * (teta * (dnorm(delta.st) * p2))))) * ((1/teta) * 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * ((1 - 
        (pnorm(delta.st) + epsilon) * p1)^(teta - 1) * (teta * 
        (dnorm(delta.st) * p1))) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 
        1) * (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p1)^teta)) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^teta) + (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^teta) * ((1 - 
        (pnorm(delta.st) + epsilon) * p2)^(teta - 1) * (teta * 
        (dnorm(delta.st) * p2)))))) - dnorm(delta.st)/(pnorm(delta.st) + 
    epsilon)^2 * ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^teta) * (1 - 
    (1 - (pnorm(delta.st) + epsilon) * p2)^teta))^((1/teta) - 
    1) * ((1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
    ((1 - (pnorm(delta.st) + epsilon) * p1)^(teta - 1) * (teta * 
        (dnorm(delta.st) * p1))) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 
    1) * (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p1)^teta)) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    p2)^teta) + (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^teta) * ((1 - 
    (pnorm(delta.st) + epsilon) * p2)^(teta - 1) * (teta * (dnorm(delta.st) * 
    p2)))))) - (dnorm(delta.st)/(pnorm(delta.st) + epsilon)^2 * 
    ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * p1)^teta) * (1 - (1 - 
        (pnorm(delta.st) + epsilon) * p2)^teta))^((1/teta) - 
        1) * ((1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
        ((1 - (pnorm(delta.st) + epsilon) * p1)^(teta - 1) * 
            (teta * (dnorm(delta.st) * p1))) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 
        1) * (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p1)^teta)) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^teta) + (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^teta) * ((1 - 
        (pnorm(delta.st) + epsilon) * p2)^(teta - 1) * (teta * 
        (dnorm(delta.st) * p2)))))) - (delta.st * dnorm(delta.st)/(pnorm(delta.st) + 
    epsilon)^2 + dnorm(delta.st) * (2 * (dnorm(delta.st) * (pnorm(delta.st) + 
    epsilon)))/((pnorm(delta.st) + epsilon)^2)^2) * (1 - (1 - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - (1 - 
        (pnorm(delta.st) + epsilon) * p1)^teta) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^teta))^(1/teta))) 




c.copula2.be1del <-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * (((1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 
        1) * teta - (1 - delta * p1)^((teta - 1) - 1) * ((teta - 
        1) * p1) * (teta * delta)) - (1 - (1 - delta)^teta)^-(1 + 
        1) * ((1 - delta)^(teta - 1) * teta) * ((1 - delta * 
        p1)^(teta - 1) * (teta * delta))) * (1 - (1 - delta * 
        p2)^teta) + (1 - (1 - delta)^teta)^-1 * ((1 - delta * 
        p1)^(teta - 1) * (teta * delta)) * ((1 - delta * p2)^(teta - 
        1) * (teta * p2)))) - (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * (1 - (1 - delta * p2)^teta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * (((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * p1)^(teta - 1) * (teta * p1)) - (1 - (1 - delta)^teta)^-(1 + 
    1) * ((1 - delta)^(teta - 1) * teta) * (1 - (1 - delta * 
    p1)^teta)) * (1 - (1 - delta * p2)^teta) + (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * ((1 - delta * p2)^(teta - 1) * 
    (teta * p2)))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * p1)^(teta - 1) * (teta * delta)) * (1 - (1 - 
    delta * p2)^teta)))) - 1/delta^2 * ((1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 
    1) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * ((1 - delta * 
    p1)^(teta - 1) * (teta * delta)) * (1 - (1 - delta * p2)^teta)))))*dnorm(delta.st)

        
        


        
c.copula2.be2del <- (1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    p1)^teta) * (1 - (1 - delta * p2)^teta))^((1/teta) - 1) * 
    ((1/teta) * (((1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 
        1) * (teta * p1)) - (1 - (1 - delta)^teta)^-(1 + 1) * 
        ((1 - delta)^(teta - 1) * teta) * (1 - (1 - delta * p1)^teta)) * 
        ((1 - delta * p2)^(teta - 1) * (teta * delta)) + (1 - 
        (1 - delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * 
        ((1 - delta * p2)^(teta - 1) * teta - (1 - delta * p2)^((teta - 
            1) - 1) * ((teta - 1) * p2) * (teta * delta)))) - 
    (1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * 
        (1 - (1 - delta * p2)^teta))^(((1/teta) - 1) - 1) * (((1/teta) - 
        1) * (((1 - (1 - delta)^teta)^-1 * ((1 - delta * p1)^(teta - 
        1) * (teta * p1)) - (1 - (1 - delta)^teta)^-(1 + 1) * 
        ((1 - delta)^(teta - 1) * teta) * (1 - (1 - delta * p1)^teta)) * 
        (1 - (1 - delta * p2)^teta) + (1 - (1 - delta)^teta)^-1 * 
        (1 - (1 - delta * p1)^teta) * ((1 - delta * p2)^(teta - 
        1) * (teta * p2)))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
        (1 - (1 - delta * p1)^teta) * ((1 - delta * p2)^(teta - 
        1) * (teta * delta))))) - 1/delta^2 * ((1 - (1 - (1 - 
    delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * (1 - (1 - 
    delta * p2)^teta))^((1/teta) - 1) * ((1/teta) * ((1 - (1 - 
    delta)^teta)^-1 * (1 - (1 - delta * p1)^teta) * ((1 - delta * 
    p2)^(teta - 1) * (teta * delta))))))*dnorm(delta.st)
 



bit1.thdel <-1/(pnorm(delta.st) + epsilon) * ((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p1)^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
    epsilon))) * exp(teta.st))) * ((1 - (pnorm(delta.st) + epsilon) * 
    p1)^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * 
    p1))) - ((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1) * (dnorm(delta.st)/(1 - (pnorm(delta.st) + epsilon)) * 
    exp(teta.st)) + (1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st)) * (log((1 - 
    (pnorm(delta.st) + epsilon))) * exp(teta.st))) + (1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 1 + 
    1) * ((1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st)))) * ((1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1) * (log((1 - 
    (pnorm(delta.st) + epsilon))) * exp(teta.st)))) * (1 - (1 - 
    (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 1)) + ((1 - 
    (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
    ((1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 1) * 
        (dnorm(delta.st) * p1/(1 - (pnorm(delta.st) + epsilon) * 
            p1) * exp(teta.st)) + (1 - (pnorm(delta.st) + epsilon) * 
        p1)^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
        (dnorm(delta.st) * p1)) * (log((1 - (pnorm(delta.st) + 
        epsilon) * p1)) * exp(teta.st))) + (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
    dnorm(delta.st))) * ((1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
    1) * (log((1 - (pnorm(delta.st) + epsilon) * p1)) * exp(teta.st))))) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1)) + ((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1) * (log((1 - (pnorm(delta.st) + epsilon))) * exp(teta.st))) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
        1)) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
    1) * (log((1 - (pnorm(delta.st) + epsilon) * p1)) * exp(teta.st)))) * 
    ((1 - (pnorm(delta.st) + epsilon) * p2)^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p2))) - 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p1))) - 
        (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
            1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
            1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * 
            (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
                1))) * ((1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon) * p2)) * exp(teta.st))) - 
        (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
            1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
            p1)^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
            epsilon) * p2)^(exp(teta.st) + 1) * (dnorm(delta.st) * 
            p2/(1 - (pnorm(delta.st) + epsilon) * p2) * exp(teta.st)) + 
            (1 - (pnorm(delta.st) + epsilon) * p2)^((exp(teta.st) + 
                1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * 
                p2)) * (log((1 - (pnorm(delta.st) + epsilon) * 
                p2)) * exp(teta.st)))))) - (1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p1)^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^(exp(teta.st) + 1)))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * ((1 - 
    (pnorm(delta.st) + epsilon) * p1)^((exp(teta.st) + 1) - 1) * 
    ((exp(teta.st) + 1) * (dnorm(delta.st) * p1))) - (1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 1) * 
    ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 1) - 
        1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - (1 - 
    (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 1))) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1)) + (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
    1)) * ((1 - (pnorm(delta.st) + epsilon) * p2)^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p2))))) * 
    ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon))) * exp(teta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p1)^(exp(teta.st) + 1)) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-1 * ((1 - (pnorm(delta.st) + 
        epsilon) * p1)^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon) * p1)) * exp(teta.st)))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1)) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p1)^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon) * p2)) * exp(teta.st))))) - ((1 - (1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * (1 - 
    (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 1)) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1)))^(1/(exp(teta.st) + 1)) * ((((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * ((1 - (pnorm(delta.st) + 
    epsilon) * p1)^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
    1) * (dnorm(delta.st) * p1))) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
    dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    p1)^(exp(teta.st) + 1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    p2)^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
    1)) * ((1 - (pnorm(delta.st) + epsilon) * p2)^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p2))))/(1 - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
            1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1))) * (exp(teta.st)/(exp(teta.st) + 1)^2)) + (1 - (1 - 
    (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
        1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
    1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p1))) - 
        (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
            1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
            1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * 
            (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
                1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
        p2)^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p1)^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
        epsilon) * p2)^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
        1) * (dnorm(delta.st) * p2))))) * (log((1 - (1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * (1 - 
    (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 1)) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1)))) * (exp(teta.st)/(exp(teta.st) + 1)^2)))) - dnorm(delta.st)/(pnorm(delta.st) + 
    epsilon)^2 * ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
    1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
    1)))^(1/(exp(teta.st) + 1)) * (log((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p1)^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
    1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
    1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon))) * exp(teta.st))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
            1)) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon) * p1)) * exp(teta.st)))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
            1)) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
        1)) * ((1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon) * p2)) * exp(teta.st))))))


}


if(BivD=="BB8.90"){

epsilon <- 0
   
  c.copula.be1 <--1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
    (1 - p1))^-(teta + 1) * (teta * delta)) * (1 - (1 + delta * 
    p2)^-teta))))


 
  c.copula.be2 <-1 - -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
    delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 1) * 
    (teta * delta)))))

  c.copula.theta <- (-1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^(-1/teta) * 
    (log((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
        (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))) * (1/teta^2)) - 
    (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * (1 - 
        p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
        1) * ((-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + 
        delta * (1 - p1))^-teta * log((1 + delta * (1 - p1)))) - 
        (1 - (1 + delta)^-teta)^-(1 + 1) * ((1 + delta)^-teta * 
            log((1 + delta))) * (1 - (1 + delta * (1 - p1))^-teta)) * 
        (1 - (1 + delta * p2)^-teta) + (1 - (1 + delta)^-teta)^-1 * 
        (1 - (1 + delta * (1 - p1))^-teta) * ((1 + delta * p2)^-teta * 
        log((1 + delta * p2)))))))*(-exp(teta.st))


  c.copula.delta <- (-(1/delta^2 * (1 - (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
    delta * (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^(-1/teta)) + 
    -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
        (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
        1) * ((-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + 
        delta * (1 - p1))^-(teta + 1) * (teta * (1 - p1))) - 
        (1 - (1 + delta)^-teta)^-(1 + 1) * ((1 + delta)^-(teta + 
            1) * teta) * (1 - (1 + delta * (1 - p1))^-teta)) * 
        (1 - (1 + delta * p2)^-teta) + (1 - (1 + delta)^-teta)^-1 * 
        (1 - (1 + delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 
        1) * (teta * p2)))))))*(-dnorm(delta.st))


 

 c.copula2.be1 <- -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((1 - (1 + delta)^-teta)^-1 * 
    ((1 + delta * (1 - p1))^-(teta + 1) * (teta * delta)) * (1 - 
    (1 + delta * p2)^-teta))) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * 
    ((1 + delta * (1 - p1))^-(teta + 1) * (teta * delta)) * (1 - 
    (1 + delta * p2)^-teta))) + (1 - (1 - (1 + delta)^-teta)^-1 * 
    (1 - (1 + delta * (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
    (1 - p1))^-(teta + 1 + 1) * ((teta + 1) * delta) * (teta * 
    delta)) * (1 - (1 + delta * p2)^-teta))))

                  
 c.copula2.be2 <--1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
    delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 1 + 
    1) * ((teta + 1) * delta) * (teta * delta)))) + (1 - (1 - 
    (1 + delta)^-teta)^-1 * (1 - (1 + delta * (1 - p1))^-teta) * 
    (1 - (1 + delta * p2)^-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * (1 - 
    p1))^-teta) * ((1 + delta * p2)^-(teta + 1) * (teta * delta)))) * 
    ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
        (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 1) * (teta * 
        delta)))))

             




c.copula2.be1be2 <- -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
    (1 - p1))^-(teta + 1) * (teta * delta)) * ((1 + delta * p2)^-(teta + 
    1) * (teta * delta)))) - (1 - (1 - (1 + delta)^-teta)^-1 * 
    (1 - (1 + delta * (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * ((1 - (1 + delta)^-teta)^-1 * 
    (1 - (1 + delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 
    1) * (teta * delta)))) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * 
    ((1 + delta * (1 - p1))^-(teta + 1) * (teta * delta)) * (1 - 
    (1 + delta * p2)^-teta))))






c.copula2.be1th <-(-1/delta * (((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * (log((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))) * (1/teta^2)) - 
    (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * (1 - 
        p1))^-teta) * (1 - (1 + delta * p2)^-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * (((1 - (1 + delta)^-teta)^-1 * 
        ((1 + delta * (1 - p1))^-teta * log((1 + delta * (1 - 
            p1)))) - (1 - (1 + delta)^-teta)^-(1 + 1) * ((1 + 
        delta)^-teta * log((1 + delta))) * (1 - (1 + delta * 
        (1 - p1))^-teta)) * (1 - (1 + delta * p2)^-teta) + (1 - 
        (1 + delta)^-teta)^-1 * (1 - (1 + delta * (1 - p1))^-teta) * 
        ((1 + delta * p2)^-teta * log((1 + delta * p2)))))) * 
    ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
        (1 - p1))^-(teta + 1) * (teta * delta)) * (1 - (1 + delta * 
        p2)^-teta))) + (1 - (1 - (1 + delta)^-teta)^-1 * (1 - 
    (1 + delta * (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * (1/teta^2 * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
    (1 - p1))^-(teta + 1) * (teta * delta)) * (1 - (1 + delta * 
    p2)^-teta)) + (-1/teta) * (((1 - (1 + delta)^-teta)^-1 * 
    ((1 + delta * (1 - p1))^-(teta + 1) * delta - (1 + delta * 
        (1 - p1))^-(teta + 1) * log((1 + delta * (1 - p1))) * 
        (teta * delta)) - (1 - (1 + delta)^-teta)^-(1 + 1) * 
    ((1 + delta)^-teta * log((1 + delta))) * ((1 + delta * (1 - 
    p1))^-(teta + 1) * (teta * delta))) * (1 - (1 + delta * p2)^-teta) + 
    (1 - (1 + delta)^-teta)^-1 * ((1 + delta * (1 - p1))^-(teta + 
        1) * (teta * delta)) * ((1 + delta * p2)^-teta * log((1 + 
        delta * p2)))))))*(-exp(teta.st))


        


c.copula2.be2th <-(-(-1/delta * (((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * (log((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))) * (1/teta^2)) - 
    (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * (1 - 
        p1))^-teta) * (1 - (1 + delta * p2)^-teta))^(((-1/teta) - 
        1) - 1) * (((-1/teta) - 1) * (((1 - (1 + delta)^-teta)^-1 * 
        ((1 + delta * (1 - p1))^-teta * log((1 + delta * (1 - 
            p1)))) - (1 - (1 + delta)^-teta)^-(1 + 1) * ((1 + 
        delta)^-teta * log((1 + delta))) * (1 - (1 + delta * 
        (1 - p1))^-teta)) * (1 - (1 + delta * p2)^-teta) + (1 - 
        (1 + delta)^-teta)^-1 * (1 - (1 + delta * (1 - p1))^-teta) * 
        ((1 + delta * p2)^-teta * log((1 + delta * p2)))))) * 
    ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
        (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 1) * (teta * 
        delta)))) + (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
    delta * (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * (1/teta^2 * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
    delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 1) * 
    (teta * delta))) + (-1/teta) * (((1 - (1 + delta)^-teta)^-1 * 
    ((1 + delta * (1 - p1))^-teta * log((1 + delta * (1 - p1)))) - 
    (1 - (1 + delta)^-teta)^-(1 + 1) * ((1 + delta)^-teta * log((1 + 
        delta))) * (1 - (1 + delta * (1 - p1))^-teta)) * ((1 + 
    delta * p2)^-(teta + 1) * (teta * delta)) + (1 - (1 + delta)^-teta)^-1 * 
    (1 - (1 + delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 
    1) * delta - (1 + delta * p2)^-(teta + 1) * log((1 + delta * 
    p2)) * (teta * delta)))))))*(-exp(teta.st))








bit1.th2 <--(-1/delta * ((1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + 
    delta * p2)^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 1)) * 
    (log((1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - 
        (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + 
        delta * p2)^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
        1)))/((exp(teta.st) + 1)^2)^2) - (((1 - (1 + delta)^(exp(teta.st) + 
        1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
        delta)) * exp(teta.st))) * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 
        1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + 
        delta * (1 - p1))^(exp(teta.st) + 1) * (log((1 + delta * 
        (1 - p1))) * exp(teta.st)))) * (1 - (1 + delta * p2)^(exp(teta.st) + 
        1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - 
        (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * ((1 + delta * 
        p2)^(exp(teta.st) + 1) * (log((1 + delta * p2)) * exp(teta.st))))/(1 - 
        (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * 
            (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + delta * 
            p2)^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)) - ((1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + 
    delta * p2)^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 1)) * 
    (log((1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - 
        (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + 
        delta * p2)^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)) + (1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + 
    delta * p2)^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 + delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
    delta)) * exp(teta.st))) * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 
    1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 + delta * (1 - p1))) * 
    exp(teta.st)))) * (1 - (1 + delta * p2)^(exp(teta.st) + 1)) - 
    (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * 
        (1 - p1))^(exp(teta.st) + 1)) * ((1 + delta * p2)^(exp(teta.st) + 
        1) * (log((1 + delta * p2)) * exp(teta.st)))))) * (log((1 - 
    (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * 
        (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + delta * p2)^(exp(teta.st) + 
        1)))) * (exp(teta.st)/(exp(teta.st) + 1)^2)) + ((1 - 
    (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * 
        (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + delta * p2)^(exp(teta.st) + 
        1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 
    1)) * ((((1 - (1 + delta)^(exp(teta.st) + 1))^-(1 + 1 + 1) * 
    ((1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + delta)) * 
        exp(teta.st)))) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
    delta)) * exp(teta.st))) + (1 - (1 + delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
    delta)) * exp(teta.st)) * (log((1 + delta)) * exp(teta.st)) + 
    (1 + delta)^(exp(teta.st) + 1) * (log((1 + delta)) * exp(teta.st)))) * 
    (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 + 
    delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 
    1) * (log((1 + delta)) * exp(teta.st))) * ((1 + delta * (1 - 
    p1))^(exp(teta.st) + 1) * (log((1 + delta * (1 - p1))) * 
    exp(teta.st))) - ((1 - (1 + delta)^(exp(teta.st) + 1))^-(1 + 
    1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + delta)) * 
    exp(teta.st))) * ((1 + delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 + delta * (1 - p1))) * exp(teta.st))) + (1 - 
    (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 + delta * (1 - p1))) * exp(teta.st)) * (log((1 + 
    delta * (1 - p1))) * exp(teta.st)) + (1 + delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 + delta * (1 - p1))) * exp(teta.st))))) * (1 - 
    (1 + delta * p2)^(exp(teta.st) + 1)) - ((1 - (1 + delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
    delta)) * exp(teta.st))) * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 
    1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 + delta * (1 - p1))) * 
    exp(teta.st)))) * ((1 + delta * p2)^(exp(teta.st) + 1) * 
    (log((1 + delta * p2)) * exp(teta.st))) - (((1 - (1 + delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
    delta)) * exp(teta.st))) * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 
    1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 + delta * (1 - p1))) * 
    exp(teta.st)))) * ((1 + delta * p2)^(exp(teta.st) + 1) * 
    (log((1 + delta * p2)) * exp(teta.st))) + (1 - (1 + delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * 
    ((1 + delta * p2)^(exp(teta.st) + 1) * (log((1 + delta * 
        p2)) * exp(teta.st)) * (log((1 + delta * p2)) * exp(teta.st)) + 
        (1 + delta * p2)^(exp(teta.st) + 1) * (log((1 + delta * 
            p2)) * exp(teta.st))))) - exp(teta.st)/(exp(teta.st) + 
    1)^2 * (((1 - (1 + delta)^(exp(teta.st) + 1))^-(1 + 1) * 
    ((1 + delta)^(exp(teta.st) + 1) * (log((1 + delta)) * exp(teta.st))) * 
    (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 + 
    delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 + delta * (1 - p1))) * exp(teta.st)))) * (1 - 
    (1 + delta * p2)^(exp(teta.st) + 1)) - (1 - (1 + delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * 
    ((1 + delta * p2)^(exp(teta.st) + 1) * (log((1 + delta * 
        p2)) * exp(teta.st))))) - ((1 - (1 - (1 + delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * 
    (1 - (1 + delta * p2)^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 
    1)) - 1) * (log((1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + 
    delta * p2)^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - 
    (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 + delta * 
    p2)^(exp(teta.st) + 1)))^(((1/(exp(teta.st) + 1)) - 1) - 
    1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 + delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
    delta)) * exp(teta.st))) * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 
    1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 + delta * (1 - p1))) * 
    exp(teta.st)))) * (1 - (1 + delta * p2)^(exp(teta.st) + 1)) - 
    (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * 
        (1 - p1))^(exp(teta.st) + 1)) * ((1 + delta * p2)^(exp(teta.st) + 
        1) * (log((1 + delta * p2)) * exp(teta.st)))))) * ((1/(exp(teta.st) + 
    1)) * (((1 - (1 + delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 + 
    delta)^(exp(teta.st) + 1) * (log((1 + delta)) * exp(teta.st))) * 
    (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 + 
    delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 + delta * (1 - p1))) * exp(teta.st)))) * (1 - 
    (1 + delta * p2)^(exp(teta.st) + 1)) - (1 - (1 + delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 + delta * (1 - p1))^(exp(teta.st) + 1)) * 
    ((1 + delta * p2)^(exp(teta.st) + 1) * (log((1 + delta * 
        p2)) * exp(teta.st))))))))





bit1.del2 <--(1/(pnorm(delta.st) + epsilon) * ((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^-teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^-teta) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    p2)^-teta))^((-1/teta) - 1) * ((-1/teta) * ((((1 - (1 - (pnorm(delta.st) + 
    epsilon))^-teta)^-(1 + 1 + 1) * ((1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^-(teta + 1) * (teta * dnorm(delta.st)))) * ((1 - 
    (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * dnorm(delta.st))) + 
    (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 1) * 
        ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1 + 1) * 
            ((teta + 1) * dnorm(delta.st)) * (teta * dnorm(delta.st)) - 
            (1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
                (delta.st * dnorm(delta.st))))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^-teta) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^-teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 
    1) * (teta * dnorm(delta.st))) * ((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^-(teta + 1) * (teta * (dnorm(delta.st) * 
    (1 - p1)))) - ((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 
    1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
    dnorm(delta.st))) * ((1 - (pnorm(delta.st) + epsilon) * (1 - 
    p1))^-(teta + 1) * (teta * (dnorm(delta.st) * (1 - p1)))) + 
    (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * ((1 - 
        (pnorm(delta.st) + epsilon) * (1 - p1))^-(teta + 1 + 
        1) * ((teta + 1) * (dnorm(delta.st) * (1 - p1))) * (teta * 
        (dnorm(delta.st) * (1 - p1))) - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^-(teta + 1) * (teta * (delta.st * 
        dnorm(delta.st) * (1 - p1)))))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^-teta) - ((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 
    1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
    dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^-teta) - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
    ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-(teta + 1) * 
        (teta * (dnorm(delta.st) * (1 - p1))))) * ((1 - (pnorm(delta.st) + 
    epsilon) * p2)^-(teta + 1) * (teta * (dnorm(delta.st) * p2))) - 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 1) * 
        ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
            dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^-teta) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^-teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
        (1 - p1))^-(teta + 1) * (teta * (dnorm(delta.st) * (1 - 
        p1))))) * ((1 - (pnorm(delta.st) + epsilon) * p2)^-(teta + 
        1) * (teta * (dnorm(delta.st) * p2))) + (1 - (1 - (pnorm(delta.st) + 
        epsilon))^-teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
        (1 - p1))^-teta) * ((1 - (pnorm(delta.st) + epsilon) * 
        p2)^-(teta + 1 + 1) * ((teta + 1) * (dnorm(delta.st) * 
        p2)) * (teta * (dnorm(delta.st) * p2)) - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^-(teta + 1) * (teta * (delta.st * dnorm(delta.st) * 
        p2)))))) - (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-teta) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^-teta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((1 - (1 - (pnorm(delta.st) + 
    epsilon))^-teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 
    1) * (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^-teta) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^-teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^-(teta + 1) * (teta * (dnorm(delta.st) * (1 - p1))))) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^-teta) - (1 - 
    (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^-teta) * ((1 - (pnorm(delta.st) + epsilon) * 
    p2)^-(teta + 1) * (teta * (dnorm(delta.st) * p2))))) * ((-1/teta) * 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 1) * 
        ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
            dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^-teta) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^-teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
        (1 - p1))^-(teta + 1) * (teta * (dnorm(delta.st) * (1 - 
        p1))))) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^-teta) - 
        (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * (1 - 
            (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-teta) * 
            ((1 - (pnorm(delta.st) + epsilon) * p2)^-(teta + 
                1) * (teta * (dnorm(delta.st) * p2)))))) - dnorm(delta.st)/(pnorm(delta.st) + 
    epsilon)^2 * ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-teta) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 
    1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
    dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^-teta) - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
    ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-(teta + 1) * 
        (teta * (dnorm(delta.st) * (1 - p1))))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^-teta) - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-teta) * 
    ((1 - (pnorm(delta.st) + epsilon) * p2)^-(teta + 1) * (teta * 
        (dnorm(delta.st) * p2)))))) - (dnorm(delta.st)/(pnorm(delta.st) + 
    epsilon)^2 * ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-teta) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 
    1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
    dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^-teta) - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
    ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-(teta + 1) * 
        (teta * (dnorm(delta.st) * (1 - p1))))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^-teta) - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^-teta) * 
    ((1 - (pnorm(delta.st) + epsilon) * p2)^-(teta + 1) * (teta * 
        (dnorm(delta.st) * p2)))))) - (delta.st * dnorm(delta.st)/(pnorm(delta.st) + 
    epsilon)^2 + dnorm(delta.st) * (2 * (dnorm(delta.st) * (pnorm(delta.st) + 
    epsilon)))/((pnorm(delta.st) + epsilon)^2)^2) * (1 - (1 - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * (1 - (1 - 
        (pnorm(delta.st) + epsilon) * (1 - p1))^-teta) * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * p2)^-teta))^(-1/teta))))




c.copula2.be1del <-(1/delta^2 * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
    (1 - p1))^-(teta + 1) * (teta * delta)) * (1 - (1 + delta * 
    p2)^-teta)))) + -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * 
    (1 - (1 + delta * (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
    (1 - p1))^-(teta + 1) * teta - (1 + delta * (1 - p1))^-(teta + 
    1 + 1) * ((teta + 1) * (1 - p1)) * (teta * delta)) - (1 - 
    (1 + delta)^-teta)^-(1 + 1) * ((1 + delta)^-(teta + 1) * 
    teta) * ((1 + delta * (1 - p1))^-(teta + 1) * (teta * delta))) * 
    (1 - (1 + delta * p2)^-teta) + (1 - (1 + delta)^-teta)^-1 * 
    ((1 + delta * (1 - p1))^-(teta + 1) * (teta * delta)) * ((1 + 
    delta * p2)^-(teta + 1) * (teta * p2)))) - (1 - (1 - (1 + 
    delta)^-teta)^-1 * (1 - (1 + delta * (1 - p1))^-teta) * (1 - 
    (1 + delta * p2)^-teta))^(((-1/teta) - 1) - 1) * (((-1/teta) - 
    1) * (((1 - (1 + delta)^-teta)^-1 * ((1 + delta * (1 - p1))^-(teta + 
    1) * (teta * (1 - p1))) - (1 - (1 + delta)^-teta)^-(1 + 1) * 
    ((1 + delta)^-(teta + 1) * teta) * (1 - (1 + delta * (1 - 
    p1))^-teta)) * (1 - (1 + delta * p2)^-teta) + (1 - (1 + delta)^-teta)^-1 * 
    (1 - (1 + delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 
    1) * (teta * p2)))) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * 
    ((1 + delta * (1 - p1))^-(teta + 1) * (teta * delta)) * (1 - 
    (1 + delta * p2)^-teta)))))*(-dnorm(delta.st))


        





        
c.copula2.be2del <- (-(1/delta^2 * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
    delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 1) * 
    (teta * delta))))) + -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * 
    (1 - (1 + delta * (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^((-1/teta) - 
    1) * ((-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
    (1 - p1))^-(teta + 1) * (teta * (1 - p1))) - (1 - (1 + delta)^-teta)^-(1 + 
    1) * ((1 + delta)^-(teta + 1) * teta) * (1 - (1 + delta * 
    (1 - p1))^-teta)) * ((1 + delta * p2)^-(teta + 1) * (teta * 
    delta)) + (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
    (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 1) * teta - 
    (1 + delta * p2)^-(teta + 1 + 1) * ((teta + 1) * p2) * (teta * 
        delta)))) - (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
    delta * (1 - p1))^-teta) * (1 - (1 + delta * p2)^-teta))^(((-1/teta) - 
    1) - 1) * (((-1/teta) - 1) * (((1 - (1 + delta)^-teta)^-1 * 
    ((1 + delta * (1 - p1))^-(teta + 1) * (teta * (1 - p1))) - 
    (1 - (1 + delta)^-teta)^-(1 + 1) * ((1 + delta)^-(teta + 
        1) * teta) * (1 - (1 + delta * (1 - p1))^-teta)) * (1 - 
    (1 + delta * p2)^-teta) + (1 - (1 + delta)^-teta)^-1 * (1 - 
    (1 + delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 
    1) * (teta * p2)))) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * 
    (1 - (1 + delta * (1 - p1))^-teta) * ((1 + delta * p2)^-(teta + 
    1) * (teta * delta)))))))*(-dnorm(delta.st))




bit1.thdel <--(1/(pnorm(delta.st) + epsilon) * ((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
    epsilon))) * exp(teta.st))) * ((1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
    (dnorm(delta.st) * (1 - p1)))) - ((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1) * (dnorm(delta.st)/(1 - (pnorm(delta.st) + 
    epsilon)) * exp(teta.st)) + (1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st)) * (log((1 - 
    (pnorm(delta.st) + epsilon))) * exp(teta.st))) + (1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 1 + 
    1) * ((1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st)))) * ((1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1) * (log((1 - 
    (pnorm(delta.st) + epsilon))) * exp(teta.st)))) * (1 - (1 - 
    (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 1)) + 
    ((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
        ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1) * (dnorm(delta.st) * (1 - p1)/(1 - (pnorm(delta.st) + 
            epsilon) * (1 - p1)) * exp(teta.st)) + (1 - (pnorm(delta.st) + 
            epsilon) * (1 - p1))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
            1) * (dnorm(delta.st) * (1 - p1))) * (log((1 - (pnorm(delta.st) + 
            epsilon) * (1 - p1))) * exp(teta.st))) + (1 - (1 - 
        (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * ((1 - 
        (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - p1))) * 
        exp(teta.st))))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    p2)^(exp(teta.st) + 1)) + ((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1) * (log((1 - (pnorm(delta.st) + epsilon))) * exp(teta.st))) * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
        1)) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - p1))) * 
    exp(teta.st)))) * ((1 - (pnorm(delta.st) + epsilon) * p2)^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p2))) - 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - 
        p1)))) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1))) * ((1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon) * p2)) * exp(teta.st))) - 
        (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
            1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
            (1 - p1))^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
            epsilon) * p2)^(exp(teta.st) + 1) * (dnorm(delta.st) * 
            p2/(1 - (pnorm(delta.st) + epsilon) * p2) * exp(teta.st)) + 
            (1 - (pnorm(delta.st) + epsilon) * p2)^((exp(teta.st) + 
                1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * 
                p2)) * (log((1 - (pnorm(delta.st) + epsilon) * 
                p2)) * exp(teta.st)))))) - (1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^(exp(teta.st) + 1)))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * ((1 - 
    (pnorm(delta.st) + epsilon) * (1 - p1))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - p1)))) - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
    1)) + (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
    1)) * ((1 - (pnorm(delta.st) + epsilon) * p2)^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p2))))) * 
    ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon))) * exp(teta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 - 
        (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
        ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - 
            p1))) * exp(teta.st)))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1)) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon) * p2)) * exp(teta.st))))) - ((1 - (1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * (1 - 
    (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
        1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
    1)))^(1/(exp(teta.st) + 1)) * ((((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * ((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
    1) * (dnorm(delta.st) * (1 - p1)))) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
    dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^(exp(teta.st) + 1))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
    epsilon) * p2)^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
    1) * (dnorm(delta.st) * p2))))/(1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
    1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
    1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - 
        p1)))) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * p2)^(exp(teta.st) + 
        1)) + (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
        p1))^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + epsilon) * 
        p2)^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
        (dnorm(delta.st) * p2))))) * (log((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * p2)^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)))) - dnorm(delta.st)/(pnorm(delta.st) + epsilon)^2 * 
    ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
        p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 
        1)) * (log((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
        p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
        1)^2)) + (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
        p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 
        1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon))) * exp(teta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 - 
        (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
        ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - 
            p1))) * exp(teta.st)))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1)) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
        epsilon) * p2)^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon) * p2)) * exp(teta.st)))))))


}


if(BivD=="BB8.180"){

  
  epsilon <- 0
   
   
  c.copula.be1 <- 1 - 1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 
    1) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * ((1 - delta * 
    (1 - p1))^(teta - 1) * (teta * delta)) * (1 - (1 - delta * 
    (1 - p2))^teta))))


 
  c.copula.be2 <- 1 - 1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 
    1) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 1) * (teta * 
    delta)))))




  c.copula.theta <-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^(1/teta) * 
    (log((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - 
        p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))) * (1/teta^2)) + 
    (1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - p1))^teta) * 
        (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 1) * ((1/teta) * 
        (((1 - (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^teta * 
            log((1 - delta))) * (1 - (1 - delta * (1 - p1))^teta) - 
            (1 - (1 - delta)^teta)^-1 * ((1 - delta * (1 - p1))^teta * 
                log((1 - delta * (1 - p1))))) * (1 - (1 - delta * 
            (1 - p2))^teta) - (1 - (1 - delta)^teta)^-1 * (1 - 
            (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - 
            p2))^teta * log((1 - delta * (1 - p2))))))))*exp(teta.st)


  c.copula.delta <- (1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 
    1) * ((1/teta) * (((1 - (1 - delta)^teta)^-1 * ((1 - delta * 
    (1 - p1))^(teta - 1) * (teta * (1 - p1))) - (1 - (1 - delta)^teta)^-(1 + 
    1) * ((1 - delta)^(teta - 1) * teta) * (1 - (1 - delta * 
    (1 - p1))^teta)) * (1 - (1 - delta * (1 - p2))^teta) + (1 - 
    (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - p1))^teta) * 
    ((1 - delta * (1 - p2))^(teta - 1) * (teta * (1 - p2)))))) - 
    1/delta^2 * (1 - (1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - 
        delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^(1/teta)))*dnorm(delta.st)




 c.copula2.be1 <- -(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * (1 - p1))^(teta - 1) * (teta * delta)) * (1 - 
    (1 - delta * (1 - p2))^teta))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * (1 - p1))^(teta - 1) * (teta * delta)) * (1 - 
    (1 - delta * (1 - p2))^teta))) + (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - 
    p2))^teta))^((1/teta) - 1) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * (1 - p1))^((teta - 1) - 1) * ((teta - 1) * 
        delta) * (teta * delta)) * (1 - (1 - delta * (1 - p2))^teta)))))


                  
 c.copula2.be2 <- -(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 
    1) * (teta * delta)))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 
    1) * (teta * delta)))) + (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - 
    p2))^teta))^((1/teta) - 1) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - p2))^((teta - 
    1) - 1) * ((teta - 1) * delta) * (teta * delta))))))
                  
                  
                 




c.copula2.be1be2 <- -(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^(((1/teta) - 
    1) - 1) * (((1/teta) - 1) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 
    1) * (teta * delta)))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * (1 - p1))^(teta - 1) * (teta * delta)) * (1 - 
    (1 - delta * (1 - p2))^teta))) - (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - 
    p2))^teta))^((1/teta) - 1) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * (1 - p1))^(teta - 1) * (teta * delta)) * ((1 - 
    delta * (1 - p2))^(teta - 1) * (teta * delta))))))






c.copula2.be1th <-(-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 
    1) * ((1/teta) * (((1 - (1 - delta)^teta)^-(1 + 1) * ((1 - 
    delta)^teta * log((1 - delta))) * ((1 - delta * (1 - p1))^(teta - 
    1) * (teta * delta)) + (1 - (1 - delta)^teta)^-1 * ((1 - 
    delta * (1 - p1))^(teta - 1) * log((1 - delta * (1 - p1))) * 
    (teta * delta) + (1 - delta * (1 - p1))^(teta - 1) * delta)) * 
    (1 - (1 - delta * (1 - p2))^teta) - (1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * (1 - p1))^(teta - 1) * (teta * delta)) * ((1 - 
    delta * (1 - p2))^teta * log((1 - delta * (1 - p2))))) - 
    1/teta^2 * ((1 - (1 - delta)^teta)^-1 * ((1 - delta * (1 - 
        p1))^(teta - 1) * (teta * delta)) * (1 - (1 - delta * 
        (1 - p2))^teta))) - ((1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - 
    p2))^teta))^((1/teta) - 1) * (log((1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - 
    p2))^teta))) * (1/teta^2)) + (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - 
    p2))^teta))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^teta * log((1 - 
    delta))) * (1 - (1 - delta * (1 - p1))^teta) - (1 - (1 - 
    delta)^teta)^-1 * ((1 - delta * (1 - p1))^teta * log((1 - 
    delta * (1 - p1))))) * (1 - (1 - delta * (1 - p2))^teta) - 
    (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - p1))^teta) * 
        ((1 - delta * (1 - p2))^teta * log((1 - delta * (1 - 
            p2))))))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * (1 - p1))^(teta - 1) * (teta * delta)) * (1 - 
    (1 - delta * (1 - p2))^teta))))))*exp(teta.st)


        
        


c.copula2.be2th <-(-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 
    1) * ((1/teta) * (((1 - (1 - delta)^teta)^-(1 + 1) * ((1 - 
    delta)^teta * log((1 - delta))) * (1 - (1 - delta * (1 - 
    p1))^teta) - (1 - (1 - delta)^teta)^-1 * ((1 - delta * (1 - 
    p1))^teta * log((1 - delta * (1 - p1))))) * ((1 - delta * 
    (1 - p2))^(teta - 1) * (teta * delta)) + (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 
    1) * log((1 - delta * (1 - p2))) * (teta * delta) + (1 - 
    delta * (1 - p2))^(teta - 1) * delta)) - 1/teta^2 * ((1 - 
    (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - p1))^teta) * 
    ((1 - delta * (1 - p2))^(teta - 1) * (teta * delta)))) - 
    ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - 
        p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 
        1) * (log((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - 
        delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))) * 
        (1/teta^2)) + (1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - 
        delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * (((1 - (1 - delta)^teta)^-(1 + 
        1) * ((1 - delta)^teta * log((1 - delta))) * (1 - (1 - 
        delta * (1 - p1))^teta) - (1 - (1 - delta)^teta)^-1 * 
        ((1 - delta * (1 - p1))^teta * log((1 - delta * (1 - 
            p1))))) * (1 - (1 - delta * (1 - p2))^teta) - (1 - 
        (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - p1))^teta) * 
        ((1 - delta * (1 - p2))^teta * log((1 - delta * (1 - 
            p2))))))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
        (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - 
        p2))^(teta - 1) * (teta * delta)))))))*exp(teta.st)






bit1.th2 <-1/delta * ((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - 
    (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - delta * 
    (1 - p2))^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 1)) * (log((1 - 
    (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - delta * 
        (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - delta * (1 - 
        p2))^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2 - exp(teta.st) * (2 * (exp(teta.st) * (exp(teta.st) + 
    1)))/((exp(teta.st) + 1)^2)^2) - (((1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) * (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p1))) * 
    exp(teta.st)))) * (1 - (1 - delta * (1 - p2))^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - 
    delta * (1 - p1))^(exp(teta.st) + 1)) * ((1 - delta * (1 - 
    p2))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p2))) * 
    exp(teta.st))))/(1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - 
    delta * (1 - p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) - ((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - 
    delta * (1 - p2))^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 
    1)) * (log((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - 
    delta * (1 - p2))^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - 
    (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - delta * 
    (1 - p2))^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 1)) - 
    1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) * (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p1))) * 
    exp(teta.st)))) * (1 - (1 - delta * (1 - p2))^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - 
    delta * (1 - p1))^(exp(teta.st) + 1)) * ((1 - delta * (1 - 
    p2))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p2))) * 
    exp(teta.st)))))) * (log((1 - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * 
    (1 - (1 - delta * (1 - p2))^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + ((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - 
    delta * (1 - p2))^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * ((((1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1 + 1) * ((1 + 1) * ((1 - delta)^(exp(teta.st) + 
    1) * (log((1 - delta)) * exp(teta.st)))) * ((1 - delta)^(exp(teta.st) + 
    1) * (log((1 - delta)) * exp(teta.st))) + (1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st)) * (log((1 - delta)) * exp(teta.st)) + 
    (1 - delta)^(exp(teta.st) + 1) * (log((1 - delta)) * exp(teta.st)))) * 
    (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 - 
    delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 
    1) * (log((1 - delta)) * exp(teta.st))) * ((1 - delta * (1 - 
    p1))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p1))) * 
    exp(teta.st))) - ((1 - (1 - delta)^(exp(teta.st) + 1))^-(1 + 
    1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - delta)) * 
    exp(teta.st))) * ((1 - delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 - delta * (1 - p1))) * exp(teta.st))) + (1 - 
    (1 - delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 - delta * (1 - p1))) * exp(teta.st)) * (log((1 - 
    delta * (1 - p1))) * exp(teta.st)) + (1 - delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 - delta * (1 - p1))) * exp(teta.st))))) * (1 - 
    (1 - delta * (1 - p2))^(exp(teta.st) + 1)) - ((1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) * (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p1))) * 
    exp(teta.st)))) * ((1 - delta * (1 - p2))^(exp(teta.st) + 
    1) * (log((1 - delta * (1 - p2))) * exp(teta.st))) - (((1 - 
    (1 - delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 
    1) * (log((1 - delta)) * exp(teta.st))) * (1 - (1 - delta * 
    (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * ((1 - delta * (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
    delta * (1 - p1))) * exp(teta.st)))) * ((1 - delta * (1 - 
    p2))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p2))) * 
    exp(teta.st))) + (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * ((1 - delta * 
    (1 - p2))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p2))) * 
    exp(teta.st)) * (log((1 - delta * (1 - p2))) * exp(teta.st)) + 
    (1 - delta * (1 - p2))^(exp(teta.st) + 1) * (log((1 - delta * 
        (1 - p2))) * exp(teta.st))))) - exp(teta.st)/(exp(teta.st) + 
    1)^2 * (((1 - (1 - delta)^(exp(teta.st) + 1))^-(1 + 1) * 
    ((1 - delta)^(exp(teta.st) + 1) * (log((1 - delta)) * exp(teta.st))) * 
    (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 - 
    delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * (1 - p1))^(exp(teta.st) + 
    1) * (log((1 - delta * (1 - p1))) * exp(teta.st)))) * (1 - 
    (1 - delta * (1 - p2))^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * 
    ((1 - delta * (1 - p2))^(exp(teta.st) + 1) * (log((1 - delta * 
        (1 - p2))) * exp(teta.st))))) - ((1 - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * 
    (1 - (1 - delta * (1 - p2))^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 
    1)) - 1) * (log((1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - 
    delta * (1 - p2))^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - 
    (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - delta * 
    (1 - p2))^(exp(teta.st) + 1)))^(((1/(exp(teta.st) + 1)) - 
    1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 - delta)^(exp(teta.st) + 
    1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 1) * (log((1 - 
    delta)) * exp(teta.st))) * (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * ((1 - delta * 
    (1 - p1))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p1))) * 
    exp(teta.st)))) * (1 - (1 - delta * (1 - p2))^(exp(teta.st) + 
    1)) - (1 - (1 - delta)^(exp(teta.st) + 1))^-1 * (1 - (1 - 
    delta * (1 - p1))^(exp(teta.st) + 1)) * ((1 - delta * (1 - 
    p2))^(exp(teta.st) + 1) * (log((1 - delta * (1 - p2))) * 
    exp(teta.st)))))) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - 
    delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 - delta)^(exp(teta.st) + 
    1) * (log((1 - delta)) * exp(teta.st))) * (1 - (1 - delta * 
    (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * ((1 - delta * (1 - p1))^(exp(teta.st) + 1) * (log((1 - 
    delta * (1 - p1))) * exp(teta.st)))) * (1 - (1 - delta * 
    (1 - p2))^(exp(teta.st) + 1)) - (1 - (1 - delta)^(exp(teta.st) + 
    1))^-1 * (1 - (1 - delta * (1 - p1))^(exp(teta.st) + 1)) * 
    ((1 - delta * (1 - p2))^(exp(teta.st) + 1) * (log((1 - delta * 
        (1 - p2))) * exp(teta.st)))))))





bit1.del2 <-1/(pnorm(delta.st) + epsilon) * ((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^teta) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p2))^teta))^((1/teta) - 1) * ((1/teta) * (((1 - (1 - 
    (pnorm(delta.st) + epsilon))^teta)^-1 * ((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(teta - 1) * (teta * (dnorm(delta.st) * 
    (1 - p1)))) - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 
    1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * 
    dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^teta)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - 
    p2))^(teta - 1) * (teta * (dnorm(delta.st) * (1 - p2)))) - 
    ((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * ((1 - 
        (pnorm(delta.st) + epsilon) * (1 - p1))^(teta - 1) * 
        (teta * (delta.st * dnorm(delta.st) * (1 - p1))) + (1 - 
        (pnorm(delta.st) + epsilon) * (1 - p1))^((teta - 1) - 
        1) * ((teta - 1) * (dnorm(delta.st) * (1 - p1))) * (teta * 
        (dnorm(delta.st) * (1 - p1)))) + (1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 
        1) * (teta * dnorm(delta.st))) * ((1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(teta - 1) * (teta * (dnorm(delta.st) * 
        (1 - p1)))) + ((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 1) * 
        (teta * dnorm(delta.st))) * ((1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(teta - 1) * (teta * (dnorm(delta.st) * 
        (1 - p1)))) - ((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 1) * 
        (teta * (delta.st * dnorm(delta.st))) + (1 - (pnorm(delta.st) + 
        epsilon))^((teta - 1) - 1) * ((teta - 1) * dnorm(delta.st)) * 
        (teta * dnorm(delta.st))) + (1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-(1 + 1 + 1) * ((1 + 1) * ((1 - (pnorm(delta.st) + 
        epsilon))^(teta - 1) * (teta * dnorm(delta.st)))) * ((1 - 
        (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * dnorm(delta.st)))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^teta))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^teta) + 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * ((1 - 
        (pnorm(delta.st) + epsilon) * (1 - p1))^(teta - 1) * 
        (teta * (dnorm(delta.st) * (1 - p1)))) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 
        1) * (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^teta)) * ((1 - (pnorm(delta.st) + 
        epsilon) * (1 - p2))^(teta - 1) * (teta * (dnorm(delta.st) * 
        (1 - p2)))) - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^teta) * 
        ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(teta - 
            1) * (teta * (delta.st * dnorm(delta.st) * (1 - p2))) + 
            (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^((teta - 
                1) - 1) * ((teta - 1) * (dnorm(delta.st) * (1 - 
                p2))) * (teta * (dnorm(delta.st) * (1 - p2))))))) - 
    (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^teta) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^teta))^(((1/teta) - 
        1) - 1) * (((1/teta) - 1) * (((1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
        (1 - p1))^(teta - 1) * (teta * (dnorm(delta.st) * (1 - 
        p1)))) - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 1) * 
        (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^teta)) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p2))^teta) + (1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
        (1 - p1))^teta) * ((1 - (pnorm(delta.st) + epsilon) * 
        (1 - p2))^(teta - 1) * (teta * (dnorm(delta.st) * (1 - 
        p2)))))) * ((1/teta) * (((1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
        (1 - p1))^(teta - 1) * (teta * (dnorm(delta.st) * (1 - 
        p1)))) - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^(teta - 1) * 
        (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^teta)) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p2))^teta) + (1 - (1 - (pnorm(delta.st) + 
        epsilon))^teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
        (1 - p1))^teta) * ((1 - (pnorm(delta.st) + epsilon) * 
        (1 - p2))^(teta - 1) * (teta * (dnorm(delta.st) * (1 - 
        p2))))))) - dnorm(delta.st)/(pnorm(delta.st) + epsilon)^2 * 
    ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^teta) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^teta))^((1/teta) - 
        1) * ((1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
        ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(teta - 
            1) * (teta * (dnorm(delta.st) * (1 - p1)))) - (1 - 
        (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 1) * ((1 - 
        (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * dnorm(delta.st))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^teta)) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^teta) + 
        (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - 
            (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^teta) * 
            ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(teta - 
                1) * (teta * (dnorm(delta.st) * (1 - p2))))))) - 
    (dnorm(delta.st)/(pnorm(delta.st) + epsilon)^2 * ((1 - (1 - 
        (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - (1 - 
        (pnorm(delta.st) + epsilon) * (1 - p1))^teta) * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^teta))^((1/teta) - 
        1) * ((1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
        ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(teta - 
            1) * (teta * (dnorm(delta.st) * (1 - p1)))) - (1 - 
        (1 - (pnorm(delta.st) + epsilon))^teta)^-(1 + 1) * ((1 - 
        (pnorm(delta.st) + epsilon))^(teta - 1) * (teta * dnorm(delta.st))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^teta)) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^teta) + 
        (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * (1 - 
            (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^teta) * 
            ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(teta - 
                1) * (teta * (dnorm(delta.st) * (1 - p2))))))) - 
        (delta.st * dnorm(delta.st)/(pnorm(delta.st) + epsilon)^2 + 
            dnorm(delta.st) * (2 * (dnorm(delta.st) * (pnorm(delta.st) + 
                epsilon)))/((pnorm(delta.st) + epsilon)^2)^2) * 
            (1 - (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^teta)^-1 * 
                (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
                  p1))^teta) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
                (1 - p2))^teta))^(1/teta)))




c.copula2.be1del <- (-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 
    1) * ((1/teta) * (((1 - (1 - delta)^teta)^-1 * ((1 - delta * 
    (1 - p1))^(teta - 1) * teta - (1 - delta * (1 - p1))^((teta - 
    1) - 1) * ((teta - 1) * (1 - p1)) * (teta * delta)) - (1 - 
    (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^(teta - 1) * teta) * 
    ((1 - delta * (1 - p1))^(teta - 1) * (teta * delta))) * (1 - 
    (1 - delta * (1 - p2))^teta) + (1 - (1 - delta)^teta)^-1 * 
    ((1 - delta * (1 - p1))^(teta - 1) * (teta * delta)) * ((1 - 
    delta * (1 - p2))^(teta - 1) * (teta * (1 - p2))))) - (1 - 
    (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - p1))^teta) * 
        (1 - (1 - delta * (1 - p2))^teta))^(((1/teta) - 1) - 
    1) * (((1/teta) - 1) * (((1 - (1 - delta)^teta)^-1 * ((1 - 
    delta * (1 - p1))^(teta - 1) * (teta * (1 - p1))) - (1 - 
    (1 - delta)^teta)^-(1 + 1) * ((1 - delta)^(teta - 1) * teta) * 
    (1 - (1 - delta * (1 - p1))^teta)) * (1 - (1 - delta * (1 - 
    p2))^teta) + (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 1) * (teta * 
    (1 - p2))))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * ((1 - 
    delta * (1 - p1))^(teta - 1) * (teta * delta)) * (1 - (1 - 
    delta * (1 - p2))^teta)))) - 1/delta^2 * ((1 - (1 - (1 - 
    delta)^teta)^-1 * (1 - (1 - delta * (1 - p1))^teta) * (1 - 
    (1 - delta * (1 - p2))^teta))^((1/teta) - 1) * ((1/teta) * 
    ((1 - (1 - delta)^teta)^-1 * ((1 - delta * (1 - p1))^(teta - 
        1) * (teta * delta)) * (1 - (1 - delta * (1 - p2))^teta))))))*dnorm(delta.st)


        



        
c.copula2.be2del <- (-(1/delta * ((1 - (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * 
    (1 - p1))^teta) * (1 - (1 - delta * (1 - p2))^teta))^((1/teta) - 
    1) * ((1/teta) * (((1 - (1 - delta)^teta)^-1 * ((1 - delta * 
    (1 - p1))^(teta - 1) * (teta * (1 - p1))) - (1 - (1 - delta)^teta)^-(1 + 
    1) * ((1 - delta)^(teta - 1) * teta) * (1 - (1 - delta * 
    (1 - p1))^teta)) * ((1 - delta * (1 - p2))^(teta - 1) * (teta * 
    delta)) + (1 - (1 - delta)^teta)^-1 * (1 - (1 - delta * (1 - 
    p1))^teta) * ((1 - delta * (1 - p2))^(teta - 1) * teta - 
    (1 - delta * (1 - p2))^((teta - 1) - 1) * ((teta - 1) * (1 - 
        p2)) * (teta * delta)))) - (1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - 
    p2))^teta))^(((1/teta) - 1) - 1) * (((1/teta) - 1) * (((1 - 
    (1 - delta)^teta)^-1 * ((1 - delta * (1 - p1))^(teta - 1) * 
    (teta * (1 - p1))) - (1 - (1 - delta)^teta)^-(1 + 1) * ((1 - 
    delta)^(teta - 1) * teta) * (1 - (1 - delta * (1 - p1))^teta)) * 
    (1 - (1 - delta * (1 - p2))^teta) + (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 
    1) * (teta * (1 - p2))))) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 
    1) * (teta * delta))))) - 1/delta^2 * ((1 - (1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * (1 - (1 - delta * (1 - 
    p2))^teta))^((1/teta) - 1) * ((1/teta) * ((1 - (1 - delta)^teta)^-1 * 
    (1 - (1 - delta * (1 - p1))^teta) * ((1 - delta * (1 - p2))^(teta - 
    1) * (teta * delta)))))))*dnorm(delta.st)
  



bit1.thdel <-1/(pnorm(delta.st) + epsilon) * ((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p2))^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 
    1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
    epsilon))) * exp(teta.st))) * ((1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
    (dnorm(delta.st) * (1 - p1)))) - ((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1) * (dnorm(delta.st)/(1 - (pnorm(delta.st) + 
    epsilon)) * exp(teta.st)) + (1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st)) * (log((1 - 
    (pnorm(delta.st) + epsilon))) * exp(teta.st))) + (1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 1 + 
    1) * ((1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st)))) * ((1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1) * (log((1 - 
    (pnorm(delta.st) + epsilon))) * exp(teta.st)))) * (1 - (1 - 
    (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 1)) + 
    ((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
        ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1) * (dnorm(delta.st) * (1 - p1)/(1 - (pnorm(delta.st) + 
            epsilon) * (1 - p1)) * exp(teta.st)) + (1 - (pnorm(delta.st) + 
            epsilon) * (1 - p1))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
            1) * (dnorm(delta.st) * (1 - p1))) * (log((1 - (pnorm(delta.st) + 
            epsilon) * (1 - p1))) * exp(teta.st))) + (1 - (1 - 
        (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * ((1 - 
        (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - p1))) * 
        exp(teta.st))))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p2))^(exp(teta.st) + 1)) + ((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
    epsilon))) * exp(teta.st))) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * ((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))) * exp(teta.st)))) * ((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p2))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
    1) * (dnorm(delta.st) * (1 - p2)))) - (((1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * ((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
    1) * (dnorm(delta.st) * (1 - p1)))) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
    epsilon))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
    dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
    (1 - p1))^(exp(teta.st) + 1))) * ((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p2))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p2))) * exp(teta.st))) - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
    epsilon) * (1 - p2))^(exp(teta.st) + 1) * (dnorm(delta.st) * 
    (1 - p2)/(1 - (pnorm(delta.st) + epsilon) * (1 - p2)) * exp(teta.st)) + 
    (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - 
        p2))) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - 
        p2))) * exp(teta.st)))))) - (1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p2))^(exp(teta.st) + 1)))^(((1/(exp(teta.st) + 
    1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 - 
    (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * ((1 - 
    (pnorm(delta.st) + epsilon) * (1 - p1))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - p1)))) - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
    p2))^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
    1)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - p2)))))) * 
    ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon))) * exp(teta.st))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(exp(teta.st) + 1)) - (1 - (1 - 
        (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
        ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - 
            p1))) * exp(teta.st)))) * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p2))^(exp(teta.st) + 1)) - (1 - (1 - 
        (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - p2))) * 
        exp(teta.st))))) - ((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p2))^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 
    1)) * ((((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - p1)))) - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 
        1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
    p2))^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
    1)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^((exp(teta.st) + 
    1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - p2)))))/(1 - 
    (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
        p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
    1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
    1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - 
        p1)))) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
        1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
        (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
        p2))^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
        epsilon) * (1 - p2))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
        1) * (dnorm(delta.st) * (1 - p2)))))) * (log((1 - (1 - 
    (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
    (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
        1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
    1)))) * (exp(teta.st)/(exp(teta.st) + 1)^2)))) - dnorm(delta.st)/(pnorm(delta.st) + 
    epsilon)^2 * ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
    1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
    1)))^(1/(exp(teta.st) + 1)) * (log((1 - (1 - (1 - (pnorm(delta.st) + 
    epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p1))^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
    epsilon) * (1 - p2))^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
    1)^2)) + (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
    1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
    1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
    1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
    (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon))) * exp(teta.st))) * 
        (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
            1)) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
        1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * (1 - p1))^(exp(teta.st) + 
        1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - p1))) * 
        exp(teta.st)))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
        (1 - p2))^(exp(teta.st) + 1)) - (1 - (1 - (pnorm(delta.st) + 
        epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
        epsilon) * (1 - p1))^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
        epsilon) * (1 - p2))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
        epsilon) * (1 - p2))) * exp(teta.st))))))


}


if(BivD=="BB8.270"){
 
      
 
epsilon <- 0
        
      
     c.copula.be1 <-1 - -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
       p1)^-(teta + 1) * (teta * delta)) * (1 - (1 + delta * (1 - 
       p2))^-teta))))
   
   
   
    
     c.copula.be2 <- -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
       delta * p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 1) * 
       (teta * delta)))))
   
   
     c.copula.theta <- (-1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^(-1/teta) * 
       (log((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
           p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))) * (1/teta^2)) - 
       (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * 
           (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 1) * 
           ((-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
               p1)^-teta * log((1 + delta * p1))) - (1 - (1 + delta)^-teta)^-(1 + 
               1) * ((1 + delta)^-teta * log((1 + delta))) * (1 - 
               (1 + delta * p1)^-teta)) * (1 - (1 + delta * (1 - 
               p2))^-teta) + (1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
               delta * p1)^-teta) * ((1 + delta * (1 - p2))^-teta * 
               log((1 + delta * (1 - p2))))))))*(-exp(teta.st))
   
   
     c.copula.delta <- (-(1/delta^2 * (1 - (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
       delta * p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^(-1/teta)) + 
       -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
           p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
           1) * ((-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + 
           delta * p1)^-(teta + 1) * (teta * p1)) - (1 - (1 + delta)^-teta)^-(1 + 
           1) * ((1 + delta)^-(teta + 1) * teta) * (1 - (1 + delta * 
           p1)^-teta)) * (1 - (1 + delta * (1 - p2))^-teta) + (1 - 
           (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * 
           ((1 + delta * (1 - p2))^-(teta + 1) * (teta * (1 - p2))))))))*(-dnorm(delta.st))
   
   
    
   
    c.copula2.be1 <- -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
       p1)^-(teta + 1 + 1) * ((teta + 1) * delta) * (teta * delta)) * 
       (1 - (1 + delta * (1 - p2))^-teta))) + (1 - (1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^(((-1/teta) - 
       1) - 1) * (((-1/teta) - 1) * ((1 - (1 + delta)^-teta)^-1 * 
       ((1 + delta * p1)^-(teta + 1) * (teta * delta)) * (1 - (1 + 
       delta * (1 - p2))^-teta))) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * 
       ((1 + delta * p1)^-(teta + 1) * (teta * delta)) * (1 - (1 + 
       delta * (1 - p2))^-teta))))
   
   
   
                     
    c.copula2.be2 <- -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^(((-1/teta) - 
       1) - 1) * (((-1/teta) - 1) * ((1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 
       1) * (teta * delta)))) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 
       1) * (teta * delta)))) + (1 - (1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
       delta * p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 1 + 
       1) * ((teta + 1) * delta) * (teta * delta)))))
   
                     
                    
   
   
   
   c.copula2.be1be2 <- -(-1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^(((-1/teta) - 
       1) - 1) * (((-1/teta) - 1) * ((1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 
       1) * (teta * delta)))) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * 
       ((1 + delta * p1)^-(teta + 1) * (teta * delta)) * (1 - (1 + 
       delta * (1 - p2))^-teta))) - (1 - (1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
       p1)^-(teta + 1) * (teta * delta)) * ((1 + delta * (1 - p2))^-(teta + 
       1) * (teta * delta))))))
   
   
   
   
   
   c.copula2.be1th <-(-(-1/delta * (((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * (log((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))) * (1/teta^2)) - 
       (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * 
           (1 - (1 + delta * (1 - p2))^-teta))^(((-1/teta) - 1) - 
           1) * (((-1/teta) - 1) * (((1 - (1 + delta)^-teta)^-1 * 
           ((1 + delta * p1)^-teta * log((1 + delta * p1))) - (1 - 
           (1 + delta)^-teta)^-(1 + 1) * ((1 + delta)^-teta * log((1 + 
           delta))) * (1 - (1 + delta * p1)^-teta)) * (1 - (1 + 
           delta * (1 - p2))^-teta) + (1 - (1 + delta)^-teta)^-1 * 
           (1 - (1 + delta * p1)^-teta) * ((1 + delta * (1 - p2))^-teta * 
           log((1 + delta * (1 - p2))))))) * ((-1/teta) * ((1 - 
       (1 + delta)^-teta)^-1 * ((1 + delta * p1)^-(teta + 1) * (teta * 
       delta)) * (1 - (1 + delta * (1 - p2))^-teta))) + (1 - (1 - 
       (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * (1 - 
       (1 + delta * (1 - p2))^-teta))^((-1/teta) - 1) * (1/teta^2 * 
       ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * p1)^-(teta + 
           1) * (teta * delta)) * (1 - (1 + delta * (1 - p2))^-teta)) + 
       (-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
           p1)^-(teta + 1) * delta - (1 + delta * p1)^-(teta + 1) * 
           log((1 + delta * p1)) * (teta * delta)) - (1 - (1 + delta)^-teta)^-(1 + 
           1) * ((1 + delta)^-teta * log((1 + delta))) * ((1 + delta * 
           p1)^-(teta + 1) * (teta * delta))) * (1 - (1 + delta * 
           (1 - p2))^-teta) + (1 - (1 + delta)^-teta)^-1 * ((1 + 
           delta * p1)^-(teta + 1) * (teta * delta)) * ((1 + delta * 
           (1 - p2))^-teta * log((1 + delta * (1 - p2)))))))))*(-exp(teta.st))
   
   
            
   
   
   c.copula2.be2th <-(-1/delta * (((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * (log((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))) * (1/teta^2)) - 
       (1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * 
           (1 - (1 + delta * (1 - p2))^-teta))^(((-1/teta) - 1) - 
           1) * (((-1/teta) - 1) * (((1 - (1 + delta)^-teta)^-1 * 
           ((1 + delta * p1)^-teta * log((1 + delta * p1))) - (1 - 
           (1 + delta)^-teta)^-(1 + 1) * ((1 + delta)^-teta * log((1 + 
           delta))) * (1 - (1 + delta * p1)^-teta)) * (1 - (1 + 
           delta * (1 - p2))^-teta) + (1 - (1 + delta)^-teta)^-1 * 
           (1 - (1 + delta * p1)^-teta) * ((1 + delta * (1 - p2))^-teta * 
           log((1 + delta * (1 - p2))))))) * ((-1/teta) * ((1 - 
       (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * ((1 + 
       delta * (1 - p2))^-(teta + 1) * (teta * delta)))) + (1 - 
       (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * 
           (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 1) * 
       (1/teta^2 * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
           p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 1) * (teta * 
           delta))) + (-1/teta) * (((1 - (1 + delta)^-teta)^-1 * 
           ((1 + delta * p1)^-teta * log((1 + delta * p1))) - (1 - 
           (1 + delta)^-teta)^-(1 + 1) * ((1 + delta)^-teta * log((1 + 
           delta))) * (1 - (1 + delta * p1)^-teta)) * ((1 + delta * 
           (1 - p2))^-(teta + 1) * (teta * delta)) + (1 - (1 + delta)^-teta)^-1 * 
           (1 - (1 + delta * p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 
           1) * delta - (1 + delta * (1 - p2))^-(teta + 1) * log((1 + 
           delta * (1 - p2))) * (teta * delta))))))*(-exp(teta.st))
   
   
   
   
   
   
   
   bit1.th2 <--(-1/delta * ((1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * 
       (1 - (1 + delta * p1)^(exp(teta.st) + 1)) * (1 - (1 + delta * 
       (1 - p2))^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 1)) * (log((1 - 
       (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * 
           p1)^(exp(teta.st) + 1)) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 
           1)))) * (exp(teta.st)/(exp(teta.st) + 1)^2 - exp(teta.st) * 
       (2 * (exp(teta.st) * (exp(teta.st) + 1)))/((exp(teta.st) + 
       1)^2)^2) - (((1 - (1 + delta)^(exp(teta.st) + 1))^-(1 + 1) * 
       ((1 + delta)^(exp(teta.st) + 1) * (log((1 + delta)) * exp(teta.st))) * 
       (1 - (1 + delta * p1)^(exp(teta.st) + 1)) - (1 - (1 + delta)^(exp(teta.st) + 
       1))^-1 * ((1 + delta * p1)^(exp(teta.st) + 1) * (log((1 + 
       delta * p1)) * exp(teta.st)))) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 
       1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + 
       delta * p1)^(exp(teta.st) + 1)) * ((1 + delta * (1 - p2))^(exp(teta.st) + 
       1) * (log((1 + delta * (1 - p2))) * exp(teta.st))))/(1 - 
       (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * 
           p1)^(exp(teta.st) + 1)) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 
           1))) * (exp(teta.st)/(exp(teta.st) + 1)^2)) - ((1 - (1 - 
       (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * p1)^(exp(teta.st) + 
       1)) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 
       1)) * (log((1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * 
       (1 - (1 + delta * p1)^(exp(teta.st) + 1)) * (1 - (1 + delta * 
       (1 - p2))^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
       1)^2)) + (1 - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - 
       (1 + delta * p1)^(exp(teta.st) + 1)) * (1 - (1 + delta * 
       (1 - p2))^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 1)) - 
       1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 + delta)^(exp(teta.st) + 
       1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
       delta)) * exp(teta.st))) * (1 - (1 + delta * p1)^(exp(teta.st) + 
       1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * 
       p1)^(exp(teta.st) + 1) * (log((1 + delta * p1)) * exp(teta.st)))) * 
       (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 1)) - (1 - (1 + 
       delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * p1)^(exp(teta.st) + 
       1)) * ((1 + delta * (1 - p2))^(exp(teta.st) + 1) * (log((1 + 
       delta * (1 - p2))) * exp(teta.st)))))) * (log((1 - (1 - (1 + 
       delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * p1)^(exp(teta.st) + 
       1)) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 1)))) * 
       (exp(teta.st)/(exp(teta.st) + 1)^2)) + ((1 - (1 - (1 + delta)^(exp(teta.st) + 
       1))^-1 * (1 - (1 + delta * p1)^(exp(teta.st) + 1)) * (1 - 
       (1 + delta * (1 - p2))^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 
       1)) - 1) * ((1/(exp(teta.st) + 1)) * ((((1 - (1 + delta)^(exp(teta.st) + 
       1))^-(1 + 1 + 1) * ((1 + 1) * ((1 + delta)^(exp(teta.st) + 
       1) * (log((1 + delta)) * exp(teta.st)))) * ((1 + delta)^(exp(teta.st) + 
       1) * (log((1 + delta)) * exp(teta.st))) + (1 - (1 + delta)^(exp(teta.st) + 
       1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
       delta)) * exp(teta.st)) * (log((1 + delta)) * exp(teta.st)) + 
       (1 + delta)^(exp(teta.st) + 1) * (log((1 + delta)) * exp(teta.st)))) * 
       (1 - (1 + delta * p1)^(exp(teta.st) + 1)) - (1 - (1 + delta)^(exp(teta.st) + 
       1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
       delta)) * exp(teta.st))) * ((1 + delta * p1)^(exp(teta.st) + 
       1) * (log((1 + delta * p1)) * exp(teta.st))) - ((1 - (1 + 
       delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 
       1) * (log((1 + delta)) * exp(teta.st))) * ((1 + delta * p1)^(exp(teta.st) + 
       1) * (log((1 + delta * p1)) * exp(teta.st))) + (1 - (1 + 
       delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * p1)^(exp(teta.st) + 
       1) * (log((1 + delta * p1)) * exp(teta.st)) * (log((1 + delta * 
       p1)) * exp(teta.st)) + (1 + delta * p1)^(exp(teta.st) + 1) * 
       (log((1 + delta * p1)) * exp(teta.st))))) * (1 - (1 + delta * 
       (1 - p2))^(exp(teta.st) + 1)) - ((1 - (1 + delta)^(exp(teta.st) + 
       1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
       delta)) * exp(teta.st))) * (1 - (1 + delta * p1)^(exp(teta.st) + 
       1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * 
       p1)^(exp(teta.st) + 1) * (log((1 + delta * p1)) * exp(teta.st)))) * 
       ((1 + delta * (1 - p2))^(exp(teta.st) + 1) * (log((1 + delta * 
           (1 - p2))) * exp(teta.st))) - (((1 - (1 + delta)^(exp(teta.st) + 
       1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 1) * (log((1 + 
       delta)) * exp(teta.st))) * (1 - (1 + delta * p1)^(exp(teta.st) + 
       1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * ((1 + delta * 
       p1)^(exp(teta.st) + 1) * (log((1 + delta * p1)) * exp(teta.st)))) * 
       ((1 + delta * (1 - p2))^(exp(teta.st) + 1) * (log((1 + delta * 
           (1 - p2))) * exp(teta.st))) + (1 - (1 + delta)^(exp(teta.st) + 
       1))^-1 * (1 - (1 + delta * p1)^(exp(teta.st) + 1)) * ((1 + 
       delta * (1 - p2))^(exp(teta.st) + 1) * (log((1 + delta * 
       (1 - p2))) * exp(teta.st)) * (log((1 + delta * (1 - p2))) * 
       exp(teta.st)) + (1 + delta * (1 - p2))^(exp(teta.st) + 1) * 
       (log((1 + delta * (1 - p2))) * exp(teta.st))))) - exp(teta.st)/(exp(teta.st) + 
       1)^2 * (((1 - (1 + delta)^(exp(teta.st) + 1))^-(1 + 1) * 
       ((1 + delta)^(exp(teta.st) + 1) * (log((1 + delta)) * exp(teta.st))) * 
       (1 - (1 + delta * p1)^(exp(teta.st) + 1)) - (1 - (1 + delta)^(exp(teta.st) + 
       1))^-1 * ((1 + delta * p1)^(exp(teta.st) + 1) * (log((1 + 
       delta * p1)) * exp(teta.st)))) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 
       1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + 
       delta * p1)^(exp(teta.st) + 1)) * ((1 + delta * (1 - p2))^(exp(teta.st) + 
       1) * (log((1 + delta * (1 - p2))) * exp(teta.st))))) - ((1 - 
       (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * 
           p1)^(exp(teta.st) + 1)) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 
           1)))^((1/(exp(teta.st) + 1)) - 1) * (log((1 - (1 - (1 + 
       delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + delta * p1)^(exp(teta.st) + 
       1)) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 1)))) * 
       (exp(teta.st)/(exp(teta.st) + 1)^2)) + (1 - (1 - (1 + delta)^(exp(teta.st) + 
       1))^-1 * (1 - (1 + delta * p1)^(exp(teta.st) + 1)) * (1 - 
       (1 + delta * (1 - p2))^(exp(teta.st) + 1)))^(((1/(exp(teta.st) + 
       1)) - 1) - 1) * (((1/(exp(teta.st) + 1)) - 1) * (((1 - (1 + 
       delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 + delta)^(exp(teta.st) + 
       1) * (log((1 + delta)) * exp(teta.st))) * (1 - (1 + delta * 
       p1)^(exp(teta.st) + 1)) - (1 - (1 + delta)^(exp(teta.st) + 
       1))^-1 * ((1 + delta * p1)^(exp(teta.st) + 1) * (log((1 + 
       delta * p1)) * exp(teta.st)))) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 
       1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + 
       delta * p1)^(exp(teta.st) + 1)) * ((1 + delta * (1 - p2))^(exp(teta.st) + 
       1) * (log((1 + delta * (1 - p2))) * exp(teta.st)))))) * ((1/(exp(teta.st) + 
       1)) * (((1 - (1 + delta)^(exp(teta.st) + 1))^-(1 + 1) * ((1 + 
       delta)^(exp(teta.st) + 1) * (log((1 + delta)) * exp(teta.st))) * 
       (1 - (1 + delta * p1)^(exp(teta.st) + 1)) - (1 - (1 + delta)^(exp(teta.st) + 
       1))^-1 * ((1 + delta * p1)^(exp(teta.st) + 1) * (log((1 + 
       delta * p1)) * exp(teta.st)))) * (1 - (1 + delta * (1 - p2))^(exp(teta.st) + 
       1)) - (1 - (1 + delta)^(exp(teta.st) + 1))^-1 * (1 - (1 + 
       delta * p1)^(exp(teta.st) + 1)) * ((1 + delta * (1 - p2))^(exp(teta.st) + 
       1) * (log((1 + delta * (1 - p2))) * exp(teta.st))))))))
   
   
   
   
   
   bit1.del2 <--(1/(pnorm(delta.st) + epsilon) * ((1 - (1 - (1 - (pnorm(delta.st) + 
       epsilon))^-teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
       p1)^-teta) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
       p2))^-teta))^((-1/teta) - 1) * ((-1/teta) * ((((1 - (1 - 
       (pnorm(delta.st) + epsilon))^-teta)^-(1 + 1 + 1) * ((1 + 
       1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
       dnorm(delta.st)))) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 
       1) * (teta * dnorm(delta.st))) + (1 - (1 - (pnorm(delta.st) + 
       epsilon))^-teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 
       1 + 1) * ((teta + 1) * dnorm(delta.st)) * (teta * dnorm(delta.st)) - 
       (1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * (delta.st * 
           dnorm(delta.st))))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
       p1)^-teta) - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 
       1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
       dnorm(delta.st))) * ((1 - (pnorm(delta.st) + epsilon) * p1)^-(teta + 
       1) * (teta * (dnorm(delta.st) * p1))) - ((1 - (1 - (pnorm(delta.st) + 
       epsilon))^-teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 
       1) * (teta * dnorm(delta.st))) * ((1 - (pnorm(delta.st) + 
       epsilon) * p1)^-(teta + 1) * (teta * (dnorm(delta.st) * p1))) + 
       (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * ((1 - 
           (pnorm(delta.st) + epsilon) * p1)^-(teta + 1 + 1) * ((teta + 
           1) * (dnorm(delta.st) * p1)) * (teta * (dnorm(delta.st) * 
           p1)) - (1 - (pnorm(delta.st) + epsilon) * p1)^-(teta + 
           1) * (teta * (delta.st * dnorm(delta.st) * p1))))) * 
       (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-teta) - 
       ((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 1) * 
           ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
               dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
           epsilon) * p1)^-teta) - (1 - (1 - (pnorm(delta.st) + 
           epsilon))^-teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
           p1)^-(teta + 1) * (teta * (dnorm(delta.st) * p1)))) * 
           ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-(teta + 
               1) * (teta * (dnorm(delta.st) * (1 - p2)))) - (((1 - 
       (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 1) * ((1 - 
       (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * dnorm(delta.st))) * 
       (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^-teta) - (1 - 
       (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * ((1 - (pnorm(delta.st) + 
       epsilon) * p1)^-(teta + 1) * (teta * (dnorm(delta.st) * p1)))) * 
       ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-(teta + 1) * 
           (teta * (dnorm(delta.st) * (1 - p2)))) + (1 - (1 - (pnorm(delta.st) + 
       epsilon))^-teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
       p1)^-teta) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-(teta + 
       1 + 1) * ((teta + 1) * (dnorm(delta.st) * (1 - p2))) * (teta * 
       (dnorm(delta.st) * (1 - p2))) - (1 - (pnorm(delta.st) + epsilon) * 
       (1 - p2))^-(teta + 1) * (teta * (delta.st * dnorm(delta.st) * 
       (1 - p2))))))) - (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
       (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^-teta) * (1 - 
       (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-teta))^(((-1/teta) - 
       1) - 1) * (((-1/teta) - 1) * (((1 - (1 - (pnorm(delta.st) + 
       epsilon))^-teta)^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 
       1) * (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
       epsilon) * p1)^-teta) - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
       ((1 - (pnorm(delta.st) + epsilon) * p1)^-(teta + 1) * (teta * 
           (dnorm(delta.st) * p1)))) * (1 - (1 - (pnorm(delta.st) + 
       epsilon) * (1 - p2))^-teta) - (1 - (1 - (pnorm(delta.st) + 
       epsilon))^-teta)^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * 
       p1)^-teta) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-(teta + 
       1) * (teta * (dnorm(delta.st) * (1 - p2)))))) * ((-1/teta) * 
       (((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 1) * 
           ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * (teta * 
               dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
           epsilon) * p1)^-teta) - (1 - (1 - (pnorm(delta.st) + 
           epsilon))^-teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
           p1)^-(teta + 1) * (teta * (dnorm(delta.st) * p1)))) * 
           (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-teta) - 
           (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * (1 - 
               (1 - (pnorm(delta.st) + epsilon) * p1)^-teta) * ((1 - 
               (pnorm(delta.st) + epsilon) * (1 - p2))^-(teta + 
               1) * (teta * (dnorm(delta.st) * (1 - p2))))))) - 
       dnorm(delta.st)/(pnorm(delta.st) + epsilon)^2 * ((1 - (1 - 
           (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * (1 - (1 - 
           (pnorm(delta.st) + epsilon) * p1)^-teta) * (1 - (1 - 
           (pnorm(delta.st) + epsilon) * (1 - p2))^-teta))^((-1/teta) - 
           1) * ((-1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 
           1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * 
           (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
           epsilon) * p1)^-teta) - (1 - (1 - (pnorm(delta.st) + 
           epsilon))^-teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
           p1)^-(teta + 1) * (teta * (dnorm(delta.st) * p1)))) * 
           (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-teta) - 
           (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * (1 - 
               (1 - (pnorm(delta.st) + epsilon) * p1)^-teta) * ((1 - 
               (pnorm(delta.st) + epsilon) * (1 - p2))^-(teta + 
               1) * (teta * (dnorm(delta.st) * (1 - p2))))))) - 
       (dnorm(delta.st)/(pnorm(delta.st) + epsilon)^2 * ((1 - (1 - 
           (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * (1 - (1 - 
           (pnorm(delta.st) + epsilon) * p1)^-teta) * (1 - (1 - 
           (pnorm(delta.st) + epsilon) * (1 - p2))^-teta))^((-1/teta) - 
           1) * ((-1/teta) * (((1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-(1 + 
           1) * ((1 - (pnorm(delta.st) + epsilon))^-(teta + 1) * 
           (teta * dnorm(delta.st))) * (1 - (1 - (pnorm(delta.st) + 
           epsilon) * p1)^-teta) - (1 - (1 - (pnorm(delta.st) + 
           epsilon))^-teta)^-1 * ((1 - (pnorm(delta.st) + epsilon) * 
           p1)^-(teta + 1) * (teta * (dnorm(delta.st) * p1)))) * 
           (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^-teta) - 
           (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * (1 - 
               (1 - (pnorm(delta.st) + epsilon) * p1)^-teta) * ((1 - 
               (pnorm(delta.st) + epsilon) * (1 - p2))^-(teta + 
               1) * (teta * (dnorm(delta.st) * (1 - p2))))))) - 
           (delta.st * dnorm(delta.st)/(pnorm(delta.st) + epsilon)^2 + 
               dnorm(delta.st) * (2 * (dnorm(delta.st) * (pnorm(delta.st) + 
                   epsilon)))/((pnorm(delta.st) + epsilon)^2)^2) * 
               (1 - (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^-teta)^-1 * 
                   (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^-teta) * 
                   (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
                     p2))^-teta))^(-1/teta))))
   
   
   
   
   c.copula2.be1del <-(-(1/delta^2 * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
       p1)^-(teta + 1) * (teta * delta)) * (1 - (1 + delta * (1 - 
       p2))^-teta)))) + -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
       p1)^-(teta + 1) * teta - (1 + delta * p1)^-(teta + 1 + 1) * 
       ((teta + 1) * p1) * (teta * delta)) - (1 - (1 + delta)^-teta)^-(1 + 
       1) * ((1 + delta)^-(teta + 1) * teta) * ((1 + delta * p1)^-(teta + 
       1) * (teta * delta))) * (1 - (1 + delta * (1 - p2))^-teta) + 
       (1 - (1 + delta)^-teta)^-1 * ((1 + delta * p1)^-(teta + 1) * 
           (teta * delta)) * ((1 + delta * (1 - p2))^-(teta + 1) * 
           (teta * (1 - p2))))) - (1 - (1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^(((-1/teta) - 
       1) - 1) * (((-1/teta) - 1) * (((1 - (1 + delta)^-teta)^-1 * 
       ((1 + delta * p1)^-(teta + 1) * (teta * p1)) - (1 - (1 + 
       delta)^-teta)^-(1 + 1) * ((1 + delta)^-(teta + 1) * teta) * 
       (1 - (1 + delta * p1)^-teta)) * (1 - (1 + delta * (1 - p2))^-teta) + 
       (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * 
           ((1 + delta * (1 - p2))^-(teta + 1) * (teta * (1 - p2))))) * 
       ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
           p1)^-(teta + 1) * (teta * delta)) * (1 - (1 + delta * 
           (1 - p2))^-teta))))))*(-dnorm(delta.st))
   
   
           
   
   
   
           
   c.copula2.be2del <- (1/delta^2 * ((1 - (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + 
       delta * p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 1) * 
       (teta * delta))))) + -1/delta * ((1 - (1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^((-1/teta) - 
       1) * ((-1/teta) * (((1 - (1 + delta)^-teta)^-1 * ((1 + delta * 
       p1)^-(teta + 1) * (teta * p1)) - (1 - (1 + delta)^-teta)^-(1 + 
       1) * ((1 + delta)^-(teta + 1) * teta) * (1 - (1 + delta * 
       p1)^-teta)) * ((1 + delta * (1 - p2))^-(teta + 1) * (teta * 
       delta)) + (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
       p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 1) * teta - 
       (1 + delta * (1 - p2))^-(teta + 1 + 1) * ((teta + 1) * (1 - 
           p2)) * (teta * delta)))) - (1 - (1 - (1 + delta)^-teta)^-1 * 
       (1 - (1 + delta * p1)^-teta) * (1 - (1 + delta * (1 - p2))^-teta))^(((-1/teta) - 
       1) - 1) * (((-1/teta) - 1) * (((1 - (1 + delta)^-teta)^-1 * 
       ((1 + delta * p1)^-(teta + 1) * (teta * p1)) - (1 - (1 + 
       delta)^-teta)^-(1 + 1) * ((1 + delta)^-(teta + 1) * teta) * 
       (1 - (1 + delta * p1)^-teta)) * (1 - (1 + delta * (1 - p2))^-teta) + 
       (1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * p1)^-teta) * 
           ((1 + delta * (1 - p2))^-(teta + 1) * (teta * (1 - p2))))) * 
       ((-1/teta) * ((1 - (1 + delta)^-teta)^-1 * (1 - (1 + delta * 
           p1)^-teta) * ((1 + delta * (1 - p2))^-(teta + 1) * (teta * 
           delta))))))*(-dnorm(delta.st))
   
   
   
   
   bit1.thdel <--(1/(pnorm(delta.st) + epsilon) * ((1 - (1 - (1 - (pnorm(delta.st) + 
       epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
       epsilon) * p1)^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
       epsilon) * (1 - p2))^(exp(teta.st) + 1)))^((1/(exp(teta.st) + 
       1)) - 1) * ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
       epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
       epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
       epsilon))) * exp(teta.st))) * ((1 - (pnorm(delta.st) + epsilon) * 
       p1)^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * 
       p1))) - ((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1) * (dnorm(delta.st)/(1 - (pnorm(delta.st) + epsilon)) * 
       exp(teta.st)) + (1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
       1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st)) * (log((1 - 
       (pnorm(delta.st) + epsilon))) * exp(teta.st))) + (1 - (1 - 
       (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 1 + 
       1) * ((1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
       1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st)))) * ((1 - 
       (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1) * (log((1 - 
       (pnorm(delta.st) + epsilon))) * exp(teta.st)))) * (1 - (1 - 
       (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 1)) + ((1 - 
       (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
       ((1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 1) * 
           (dnorm(delta.st) * p1/(1 - (pnorm(delta.st) + epsilon) * 
               p1) * exp(teta.st)) + (1 - (pnorm(delta.st) + epsilon) * 
           p1)^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
           (dnorm(delta.st) * p1)) * (log((1 - (pnorm(delta.st) + 
           epsilon) * p1)) * exp(teta.st))) + (1 - (1 - (pnorm(delta.st) + 
       epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
       epsilon))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 1) * 
       dnorm(delta.st))) * ((1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
       1) * (log((1 - (pnorm(delta.st) + epsilon) * p1)) * exp(teta.st))))) * 
       (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
           1)) + ((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1) * (log((1 - (pnorm(delta.st) + epsilon))) * exp(teta.st))) * 
       (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
           1)) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
       1) * (log((1 - (pnorm(delta.st) + epsilon) * p1)) * exp(teta.st)))) * 
       ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^((exp(teta.st) + 
           1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - 
           p2)))) - (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^((exp(teta.st) + 
       1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p1))) - 
       (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 
           1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
           1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
           (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
               1))) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
       1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - p2))) * 
       exp(teta.st))) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
       1)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
       1) * (dnorm(delta.st) * (1 - p2)/(1 - (pnorm(delta.st) + 
       epsilon) * (1 - p2)) * exp(teta.st)) + (1 - (pnorm(delta.st) + 
       epsilon) * (1 - p2))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
       1) * (dnorm(delta.st) * (1 - p2))) * (log((1 - (pnorm(delta.st) + 
       epsilon) * (1 - p2))) * exp(teta.st)))))) - (1 - (1 - (1 - 
       (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * (1 - 
       (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 1)) * 
       (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
           1)))^(((1/(exp(teta.st) + 1)) - 1) - 1) * (((1/(exp(teta.st) + 
       1)) - 1) * (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^((exp(teta.st) + 
       1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p1))) - 
       (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 
           1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
           1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
           (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
               1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
       p2))^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
       1)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^((exp(teta.st) + 
       1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - p2)))))) * 
       ((1/(exp(teta.st) + 1)) * (((1 - (1 - (pnorm(delta.st) + 
           epsilon))^(exp(teta.st) + 1))^-(1 + 1) * ((1 - (pnorm(delta.st) + 
           epsilon))^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
           epsilon))) * exp(teta.st))) * (1 - (1 - (pnorm(delta.st) + 
           epsilon) * p1)^(exp(teta.st) + 1)) - (1 - (1 - (pnorm(delta.st) + 
           epsilon))^(exp(teta.st) + 1))^-1 * ((1 - (pnorm(delta.st) + 
           epsilon) * p1)^(exp(teta.st) + 1) * (log((1 - (pnorm(delta.st) + 
           epsilon) * p1)) * exp(teta.st)))) * (1 - (1 - (pnorm(delta.st) + 
           epsilon) * (1 - p2))^(exp(teta.st) + 1)) - (1 - (1 - 
           (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
           (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
               1)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
           1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - p2))) * 
           exp(teta.st))))) - ((1 - (1 - (1 - (pnorm(delta.st) + 
       epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
       epsilon) * p1)^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
       epsilon) * (1 - p2))^(exp(teta.st) + 1)))^(1/(exp(teta.st) + 
       1)) * ((((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^((exp(teta.st) + 
       1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p1))) - 
       (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-(1 + 
           1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
           1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * (1 - 
           (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
               1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
       p2))^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
       1)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^((exp(teta.st) + 
       1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * (1 - p2)))))/(1 - 
       (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
           (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
               1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - 
           p2))^(exp(teta.st) + 1))) * (exp(teta.st)/(exp(teta.st) + 
       1)^2)) + (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
       1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
       1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
       (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
           1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^((exp(teta.st) + 
           1) - 1) * ((exp(teta.st) + 1) * (dnorm(delta.st) * p1))) - 
           (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
               1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^((exp(teta.st) + 
               1) - 1) * ((exp(teta.st) + 1) * dnorm(delta.st))) * 
               (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
                   1))) * (1 - (1 - (pnorm(delta.st) + epsilon) * 
           (1 - p2))^(exp(teta.st) + 1)) + (1 - (1 - (pnorm(delta.st) + 
           epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
           epsilon) * p1)^(exp(teta.st) + 1)) * ((1 - (pnorm(delta.st) + 
           epsilon) * (1 - p2))^((exp(teta.st) + 1) - 1) * ((exp(teta.st) + 
           1) * (dnorm(delta.st) * (1 - p2)))))) * (log((1 - (1 - 
       (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 1))^-1 * 
       (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
           1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
       1)))) * (exp(teta.st)/(exp(teta.st) + 1)^2)))) - dnorm(delta.st)/(pnorm(delta.st) + 
       epsilon)^2 * ((1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
       1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
       1)))^(1/(exp(teta.st) + 1)) * (log((1 - (1 - (1 - (pnorm(delta.st) + 
       epsilon))^(exp(teta.st) + 1))^-1 * (1 - (1 - (pnorm(delta.st) + 
       epsilon) * p1)^(exp(teta.st) + 1)) * (1 - (1 - (pnorm(delta.st) + 
       epsilon) * (1 - p2))^(exp(teta.st) + 1)))) * (exp(teta.st)/(exp(teta.st) + 
       1)^2)) + (1 - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
       1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
       1)) * (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
       1)))^((1/(exp(teta.st) + 1)) - 1) * ((1/(exp(teta.st) + 1)) * 
       (((1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
           1))^-(1 + 1) * ((1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
           1) * (log((1 - (pnorm(delta.st) + epsilon))) * exp(teta.st))) * 
           (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
               1)) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
           1))^-1 * ((1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
           1) * (log((1 - (pnorm(delta.st) + epsilon) * p1)) * exp(teta.st)))) * 
           (1 - (1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
               1)) - (1 - (1 - (pnorm(delta.st) + epsilon))^(exp(teta.st) + 
           1))^-1 * (1 - (1 - (pnorm(delta.st) + epsilon) * p1)^(exp(teta.st) + 
           1)) * ((1 - (pnorm(delta.st) + epsilon) * (1 - p2))^(exp(teta.st) + 
           1) * (log((1 - (pnorm(delta.st) + epsilon) * (1 - p2))) * 
           exp(teta.st)))))))



}



         list(c.copula.be1=c.copula.be1, c.copula.be2=c.copula.be2, c.copula.theta=c.copula.theta, c.copula.delta=c.copula.delta,
              c.copula2.be1=c.copula2.be1, c.copula2.be2=c.copula2.be2, c.copula2.be1be2=c.copula2.be1be2,c.copula2.be1th=c.copula2.be1th,
              c.copula2.be2th=c.copula2.be2th,bit1.th2=bit1.th2,bit1.del2=bit1.del2,c.copula2.be1del=c.copula2.be1del,
              c.copula2.be2del=c.copula2.be2del,bit1.thdel=bit1.thdel)     


}




     
























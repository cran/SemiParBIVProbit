distrHsAT <- function(y2, eta2, sigma2, margin2){

if(margin2 == "N"){

  pdf2          <- dnorm(y2, mean=eta2, sd = sqrt(sigma2))
    p2          <- pnorm(y2, mean=eta2, sd = sqrt(sigma2))
     
}


if(margin2 == "LN"){

  pdf2          <- dlnorm(y2, meanlog=eta2, sdlog = sqrt(sigma2))
    p2          <- plnorm(y2, meanlog=eta2, sdlog = sqrt(sigma2))

}


if(margin2 == "WEI"){

  pdf2          <- sqrt(sigma2)/exp(eta2)*(y2/exp(eta2))^(sqrt(sigma2)-1) * exp(-(y2/exp(eta2))^sqrt(sigma2))
                   
    p2          <-  1-exp(-(y2/exp(eta2))^sqrt(sigma2)) #  pWEI(y2,exp(eta2), sqrt(sigma2))

}


if(margin2 == "iG"){

  pdf2          <- exp(-0.5 * log(2 * pi) - log(sqrt(sigma2)) - (3/2) * log(y2) - 
                   ((y2 - exp(eta2))^2)/(2 * sigma2 * (exp(eta2)^2) * y2))
                   
    p2          <-  pnorm(((y2/exp(eta2)) - 1)/(sqrt(sigma2) * sqrt(y2))) + 
                    exp(2/(exp(eta2)*sigma2))* pnorm(-((y2/exp(eta2)) + 1)/(sqrt(sigma2) * sqrt(y2)))
                
}



if(margin2 == "LO"){

  pdf2          <- dlogis(y2,eta2,sqrt(sigma2)) # exp(-(y2-eta2)/sqrt(sigma2))/(sqrt(sigma2)*(1+exp(-(y2-eta2)/sqrt(sigma2)))^2)
    p2          <- plogis(y2,eta2,sqrt(sigma2)) #1/(1+exp(-(y2-eta2)/sqrt(sigma2)))
                
}


if(margin2 == "rGU"){

  pdf2          <- 1/sqrt(sigma2)*exp(-((y2-eta2)/sqrt(sigma2)+exp(-((y2-eta2)/sqrt(sigma2)))))
    p2          <- exp(-(exp(-(y2-eta2)/sqrt(sigma2))))
                
}



if(margin2 == "GU"){

  pdf2          <- exp(-exp((y2 - eta2)/sqrt(sigma2))) * (exp((y2 - eta2)/sqrt(sigma2)) * 
                   (1/sqrt(sigma2)))
    p2          <- 1 - exp(-exp((y2 - eta2)/sqrt(sigma2)))
    

}


if(margin2 == "GA"){



 pdf2          <-  dgamma(y2, shape = 1/sigma2, scale = exp(eta2) * sigma2)
                   
    p2          <-  pgamma(y2, shape = 1/sigma2, scale = exp(eta2) * sigma2)
    
    }
    
    
    

if(margin2 == "iGA"){


pdf2          <-  exp(1/sigma2 * eta2 + 1/sigma2 * log(1/sigma2 + 1) - lgamma(1/sigma2) - 
                     (1/sigma2 + 1) * log(y2) - ((exp(eta2) * (1/sigma2 + 1))/y2))
              
    p2          <-  1-pgamma(((exp(eta2) * (1/sigma2 + 1))/y2), shape = 1/sigma2, scale=1)
      
      
 }     















epsilon <- 0.0000001 
max.p   <- 0.9999999

  #pdf2 <- ifelse(pdf2 < epsilon, epsilon, pdf2 )

  p2 <- ifelse(p2 > max.p, max.p, p2) 
  p2 <- ifelse(p2 < epsilon, epsilon, p2) 


list(pdf2 = pdf2, p2 = p2)     


}




     
























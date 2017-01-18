pscr0 <- function(x, type = "copR"){

if(type == "copR"){

  if(x$margins[1]%in%c("N","N2"))         cat("\nMARGIN 1: Gaussian")  
  if(x$margins[1]=="GU")                  cat("\nMARGIN 1: Gumbel")    
  if(x$margins[1]=="rGU")                 cat("\nMARGIN 1: reverse Gumbel")  
  if(x$margins[1]=="LO")                  cat("\nMARGIN 1: logistic")   
  if(x$margins[1]=="LN")                  cat("\nMARGIN 1: log-normal") 
  if(x$margins[1]=="WEI")                 cat("\nMARGIN 1: Weibull") 
  if(x$margins[1]=="iG")                  cat("\nMARGIN 1: inverse Gaussian") 
  if(x$margins[1]%in%c("GA","GAi"))       cat("\nMARGIN 1: gamma")  
  if(x$margins[1]=="BE")                  cat("\nMARGIN 1: beta")    
  if(x$margins[1]=="DAGUM")               cat("\nMARGIN 1: Dagum")  
  if(x$margins[1]=="SM")                  cat("\nMARGIN 1: Singh-Maddala") 
  if(x$margins[1]=="FISK")                cat("\nMARGIN 1: Fisk") 
  if(x$margins[1]%in%c("NBI","NBIa"))     cat("\nMARGIN 1: Negative Binomial - Type I") 
  if(x$margins[1]%in%c("NBII","NBIIa"))   cat("\nMARGIN 1: Negative Binomial - Type II")  
  if(x$margins[1]=="PIG")                 cat("\nMARGIN 1: Poisson inverse Gaussian")
  if(x$margins[1]=="PO")                  cat("\nMARGIN 1: Poisson")   
  if(x$margins[1]=="ZTP")                 cat("\nMARGIN 1: Zero Truncated Poisson")    
  
  
  if(x$margins[2]%in%c("N","N2"))         cat("\nMARGIN 2: Gaussian")  
  if(x$margins[2]=="GU")     		  cat("\nMARGIN 2: Gumbel")    
  if(x$margins[2]=="rGU")    		  cat("\nMARGIN 2: reverse Gumbel")  
  if(x$margins[2]=="LO")     		  cat("\nMARGIN 2: logistic")   
  if(x$margins[2]=="LN")     		  cat("\nMARGIN 2: log-normal") 
  if(x$margins[2]=="WEI")    		  cat("\nMARGIN 2: Weibull") 
  if(x$margins[2]=="iG")     		  cat("\nMARGIN 2: inverse Gaussian") 
  if(x$margins[2]%in%c("GA","GAi")) 	  cat("\nMARGIN 2: gamma")   
  if(x$margins[2]=="BE")     		  cat("\nMARGIN 2: beta")    
  if(x$margins[2]=="DAGUM")  		  cat("\nMARGIN 2: Dagum")
  if(x$margins[2]=="SM")     		  cat("\nMARGIN 2: Singh-Maddala") 
  if(x$margins[2]=="FISK")   		  cat("\nMARGIN 2: Fisk") 
  if(x$margins[2]%in%c("NBI","NBIa"))     cat("\nMARGIN 2: Negative Binomial - Type I") 
  if(x$margins[2]%in%c("NBII","NBIIa"))   cat("\nMARGIN 2: Negative Binomial - Type II")  
  if(x$margins[2]=="PIG")                 cat("\nMARGIN 2: Poisson inverse Gaussian")  
  if(x$margins[2]=="PO")                  cat("\nMARGIN 2: Poisson")   
  if(x$margins[2]=="ZTP")                 cat("\nMARGIN 2: Zero Truncated Poisson")    
    
}


if(type == "gamls"){


  if(x$margins[1]%in%c("N","N2"))      cat("\nDistribution: Gaussian")  
  if(x$margins[1]=="GU")               cat("\nDistribution: Gumbel")    
  if(x$margins[1]=="rGU")              cat("\nDistribution: reverse Gumbel")  
  if(x$margins[1]=="LO")               cat("\nDistribution: logistic")   
  if(x$margins[1]=="LN")               cat("\nDistribution: log-normal") 
  if(x$margins[1]=="WEI")              cat("\nDistribution: Weibull") 
  if(x$margins[1]=="iG")               cat("\nDistribution: inverse Gaussian") 
  if(x$margins[1]%in%c("GA","GAi"))    cat("\nDistribution: gamma")  
  if(x$margins[1]=="FISK")             cat("\nDistribution: Fisk") 
  if(x$margins[1]=="BE")               cat("\nDistribution: beta")    
  if(x$margins[1]=="DAGUM")            cat("\nDistribution: Dagum")  
  if(x$margins[1]=="SM")               cat("\nDistribution: Singh-Maddala") 
  if(x$margins[1]%in%c("NBI","NBIa"))  cat("\nDistribution: Negative Binomial - Type I") 
  if(x$margins[1]%in%c("NBII","NBIIa"))cat("\nDistribution: Negative Binomial - Type II")
  if(x$margins[1]=="PIG")              cat("\nDistribution: Poisson inverse Gaussian") 
  if(x$margins[1]=="PO")               cat("\nDistribution: Poisson")   
  if(x$margins[1]=="ZTP")              cat("\nDistribution: Zero Truncated Poisson")    
  

}



if(type == "copSS"){

if(x$margins[1] %in% x$bl) cat("\nMARGIN 1: Bernoulli") 
if(x$margins[2] %in% x$bl) cat("\nMARGIN 2: Bernoulli") 

if(x$margins[2]%in%c("N","N2"))        cat("\nMARGIN 2: Gaussian")  
if(x$margins[2]=="GU")                 cat("\nMARGIN 2: Gumbel")    
if(x$margins[2]=="rGU")                cat("\nMARGIN 2: reverse Gumbel")  
if(x$margins[2]=="LO")                 cat("\nMARGIN 2: logistic")   
if(x$margins[2]=="LN")                 cat("\nMARGIN 2: log-normal") 
if(x$margins[2]=="WEI")                cat("\nMARGIN 2: Weibull") 
if(x$margins[2]=="iG")                 cat("\nMARGIN 2: inverse Gaussian") 
if(x$margins[2]%in%c("GA","GAi"))      cat("\nMARGIN 2: gamma")   
if(x$margins[2]=="BE")                 cat("\nMARGIN 2: beta")    
if(x$margins[2]=="DAGUM")              cat("\nMARGIN 2: Dagum")
if(x$margins[2]=="SM")                 cat("\nMARGIN 2: Singh-Maddala") 
if(x$margins[2]=="FISK")               cat("\nMARGIN 2: Fisk") 
if(x$margins[2] %in% c("NBI","NBIa"))  cat("\nMARGIN 2: Negative Binomial - Type I") 
if(x$margins[2]%in% c("NBII","NBIIa")) cat("\nMARGIN 2: Negative Binomial - Type II")
if(x$margins[2]=="PIG")    	       cat("\nMARGIN 2: Poisson inverse Gaussian") 
if(x$margins[2]=="PO")     	       cat("\nMARGIN 2: Poisson")   
if(x$margins[2]=="ZTP")    	       cat("\nMARGIN 2: Zero Truncated Poisson")  

}






         
}


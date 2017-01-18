print.gamlss <- function(x, ...){


 ppR <- pp(x)  
 
 cont1par <- ppR$cont1par
 cont2par <- ppR$cont2par
 cont3par <- ppR$cont3par
 m1l      <- ppR$m1l

  s1 <- "sigma2 = "; s1.p <- x$sigma2.a
  n1 <- "  nu = "; n1.p <- x$nu.a
  

  pscr0(x, type = "gamls")  

  cat("\n\nEQUATION 1")
  cat("\nLink function for mu:",m1l,"\n")
  cat("Formula: "); print(x$gam1$formula) 
  
  
  
  if( !is.null(x$X2) ){##
  

  if( x$margins[1] %in% cont2par){
  
  cat("\nEQUATION 2")
  if(x$margins[1] !="BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[2]])  
  
     }
  
 
if( x$margins[1] %in% cont3par){
     
  cat("\nEQUATION 2")
  cat("\nLink function for sigma2:","log","\n") 
  cat("Formula: "); print(x$formula[[2]]) 
  
  cat("\nEQUATION 3")
  cat("\nLink function for nu:","log","\n") 
  cat("Formula: "); print(x$formula[[3]]) 
     
     }
    
  }##
  
   
  
  cat("\n")
                                                                          
  if(x$margins[1] %in% cont1par) cat("n = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")  
                                                                     

  if(x$margins[1] %in% cont2par) cat( s1, format(s1.p, digits=3),
                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")   

  if(x$margins[1] %in% cont3par) cat( s1, format(s1.p, digits=3),
                                     n1, format(n1.p, digits=3), 
                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")   

invisible(x)

}


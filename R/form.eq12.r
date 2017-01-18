form.eq12 <- function(formula.eq1, data, v1, margins, m1d, m2d, copSS = FALSE, inde = NULL){
  
    y1m <- NULL
  
    formula.eq1r <- formula.eq1   
    y1 <- y1.test <- data[, v1[1]]
    
    if(copSS == TRUE) y1 <- y1.test <- y1[inde]  
       
    if( v1[1] != as.character(formula.eq1r[2]) ) y1.test <- try(data[, as.character(formula.eq1r[2])], silent = TRUE)
    if(class(y1.test) == "try-error") stop("Please check the syntax of the equations' responses.") 

    if(margins %in% c(m1d,m2d) && min(y1.test, na.rm = TRUE) < 0) stop("The response of one or both margins must be positive.")
    if(margins %in% c(m1d,m2d)){
    
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if(sum(as.numeric(is.wholenumber(y1.test))) != length(y1.test)) stop("The response of one or both margins must be discrete.")     
    }
    if(margins %in% c("ZTP") && min(y1.test, na.rm = TRUE) < 1) stop("The response of one or both margins must be greater than 0.") 
    
    if(margins %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK") && min(y1.test, na.rm = TRUE) <= 0) stop("The response of one or both margins must be positive.")
    if(margins %in% c("BE") && (min(y1.test, na.rm = TRUE) <= 0 || max(y1.test, na.rm = TRUE) >= 1) ) stop("The response of one or both margins must be in the interval (0,1).")
     
     
    # matrix useful for fitting
    
    if(margins %in% c("NBIa","NBIIa","NBI","PO","ZTP")){ # no PIG, NBII as these are numerical # a - all analytical
                                                         # default for NBII all numerical (no need for y2m)
                                                         # default for NBI half num half analyt (need y2m)
     
    ly1 <- length(y1)
    y1m <- list()
    my1 <- max(y1)
    for(i in 1:ly1){ y1m[[i]] <- seq(0, y1[i]); length(y1m[[i]]) <- my1+1} 
    y1m <- do.call(rbind, y1m)   
    
  
    if(max(y1) > 170 && margins %in% c("PO","ZTP") ) y1m <- mpfr( y1m, pmax(53, getPrec(y1))) 
    
    }
     
    if( margins %in% c("N","N2","LO","GU","rGU","GAi") ) formula.eq1 <- update(formula.eq1, (. + mean(.))/2 ~ . ) 
    if( margins %in% c(m1d, m2d) )                       formula.eq1 <- update(formula.eq1, log((. + mean(.))/2) ~ . )  
    if( margins %in% c("LN") )                           formula.eq1 <- update(formula.eq1, (log(.) + mean(log(.)))/2 ~ . )  
    if( margins %in% c("iG","GA","DAGUM","SM","FISK") )  formula.eq1 <- update(formula.eq1, log((. + mean(.))/2) ~ . )    
    if( margins %in% c("WEI") )                          formula.eq1 <- update(formula.eq1, log( exp(log(.) + 0.5772/(1.283/sqrt(var(log(.)))))  ) ~ . )     
    if( margins %in% c("BE") )                           formula.eq1 <- update(formula.eq1, qlogis((. + mean(.))/2) ~ . )    
  
  list(formula.eq1 = formula.eq1, formula.eq1r = formula.eq1r, y1 = y1, y1.test = y1.test, y1m = y1m)
  
  }
  
  
  

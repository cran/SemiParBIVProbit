startsn <- function(margins, y1){
       
log.nu.1 <- NULL      

           
        if( !(margins %in% c("GO")) ) par.est <- try( resp.check(y1, margin = margins, plots = FALSE, print.par = TRUE, i.f = TRUE), silent = TRUE)
        if( margins %in% c("GO") )    log.sig2.1 <- log((1/mean(y1))^2)              

        
        if( !(margins %in% c("GO")) ){
        
        if(class(par.est)=="try-error") {
        
 		if( margins %in% c("NBI","NBIa","PIG") )  log.sig2.1 <- log( max( (var(y1) - mean(y1))/mean(y1)^2, 0.1) )
 		if( margins %in% c("NBII","NBIIa") )      log.sig2.1 <- log( max( (var(y1)/mean(y1)) - 1, 0.1) ) 		
 		if( margins %in% c("N","LN") )            log.sig2.1 <- log(var(y1))  
 	        if( margins %in% c("N2") )                log.sig2.1 <- log(sqrt(var(y1)))  
		if( margins %in% c("LO") )                log.sig2.1 <- log( 3*var(y1)/pi^2 )   
		if( margins %in% c("iG") )                log.sig2.1 <- log( var(y1)/mean(y1)^3 )      
		if( margins %in% c("GU","rGU") )          log.sig2.1 <- log(6*var(y1)/pi^2)   
		if( margins %in% c("WEI") )               log.sig2.1 <- log( (1.283/sqrt(var(log(y1))))^2 )    
	        #if( margins %in% c("GO") )                log.sig2.1 <- log((1/mean(y1))^2)              
		if( margins %in% c("GA","GAi") )          log.sig2.1 <- log(var(y1)/mean(y1)^2)
		#if( margins %in% c("GA2") )               log.sig2.1 <- log(mean(y1)/var(y1))           
        	if( margins %in% c("DAGUM","SM","FISK") ) log.sig2.1 <- log(sqrt(2))  # log(0.01) #  log(sqrt(2))       # 0.1  
        	if( margins %in% c("BE"))                 log.sig2.1 <- qlogis( var(y1)/( mean(y1)*(1-mean(y1)) )  )        
        	
                                        } else log.sig2.1 <- par.est[2]
        

              if( margins %in% c("DAGUM","SM") ){
        	if(class(par.est)=="try-error") log.nu.1 <- log(1) else log.nu.1 <- par.est[3]
                                                }                                                                                    

        }
        
list(log.sig2.1 = log.sig2.1, log.nu.1 = log.nu.1)        


}


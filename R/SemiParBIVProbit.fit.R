SemiParBIVProbit.fit  <- function(func.opt, start.v, 
                                   rinit, rmax, iterlim, iterlimsp, pr.tolsp,
                                   respvec, VC,
                                   sp=NULL, qu.mag=NULL, naive){ 

  parsc <- rep(VC$parscale, length(start.v) )  

  fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, parscale=parsc,
                respvec=respvec, VC=VC, sp=sp, qu.mag=qu.mag, blather=TRUE, iterlim = iterlim)  

  iter.if <- fit$iterations  

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL
  
  
  if(naive == TRUE) l.sp1 <- 0 else l.sp1 <- VC$l.sp1
  


    if((VC$fp==FALSE && (l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0 || VC$l.sp6!=0)) ){

       stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0  
       

	  while( stoprule.SP > pr.tolsp ){ 


             spo <- sp 
             wor.c <- working.comp(fit) 
  
                	bs.mgfit <- try(magic(y = wor.c$Z,  
                	                      X = wor.c$X,
                	                      sp= sp, S = qu.mag$Ss,
                        	              off = qu.mag$off, rank = qu.mag$rank,
                                	      gcv = FALSE,
                                	      gamma = VC$infl.fac), silent = TRUE)
                		if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                		
                	sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1; names(sp) <- names(spo) 
                        o.ests <- c(fit$argument) 

             fit <- trust(func.opt, o.ests, rinit=rinit, rmax=rmax,  parscale=parsc,  
                          respvec=respvec, VC=VC, 
                          sp=sp, qu.mag=qu.mag, 
                          blather=TRUE, iterlim=iterlim)
                                                            
             iter.inner <- iter.inner + fit$iterations   	                                    
              
             if(iter.sp >= iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- max(abs(o.ests - c(fit$argument)))
                  
           
          }
          

       if(VC$gc.l == TRUE) gc()   
       
       #wor.c <- working.comp(fit) 
       #bs.mgfit <- magic(y = wor.c$Z,X = wor.c$X,sp= sp, S = qu.mag$Ss,
       #                  off = qu.mag$off, rank = qu.mag$rank,
       #                  gcv = FALSE, gamma = VC$infl.fac)
       
       magpp <- magic.post.proc(wor.c$X, bs.mgfit)
       

    }else{
    
    wor.c <- working.comp(fit) 
    
    bs.mgfit <- magic(wor.c$Z, wor.c$X, numeric(0), list(), numeric(0))    
    magpp    <- magic.post.proc(wor.c$X, bs.mgfit)
    
    
    
    }
    
    
    
    
    
    
    
    


                  list(fit = fit, 
                       iter.if = iter.if, 
                       conv.sp = conv.sp, 
                       iter.sp = iter.sp, 
                       iter.inner = iter.inner, 
                       bs.mgfit = bs.mgfit, 
                       wor.c = wor.c, 
                       sp = sp, magpp = magpp)
      
 
}





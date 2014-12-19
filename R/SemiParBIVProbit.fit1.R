SemiParBIVProbit.fit1  <- function(func.opt, start.v, 
                                 rinit, rmax, iterlim, iterlimsp, pr.tolsp,
                                 PL, eqPL, valPL, fitPL, 
                                 respvec, VC,
                                 sp=NULL, qu.mag=NULL){ 

spE <- sp

if( PL!="P" && (VC$l.sp1!=0 || VC$l.sp2!=0) ){
    lspe <- length(sp)
    if(eqPL=="both") exclu <- c(lspe-1,lspe) else exclu <- c(lspe)  
    spE <- sp[-exclu]
    qu.mag$exclu <- exclu
}

  fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax,
                sp.xi1=sp["xi1"], sp.xi2=sp["xi2"], PL=PL, eqPL=eqPL, 
                valPL=valPL, fitPL=fitPL, respvec=respvec, VC=VC,
                sp=spE, qu.mag=qu.mag, 
                blather=TRUE, iterlim = iterlim)  

  iter.if <- fit$iterations  

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL

if( PL!="P" && fitPL!="pLiksp" ) list(fit=fit, iter.if=iter.if, conv.sp=conv.sp, iter.sp=iter.sp, iter.inner=iter.inner, bs.mgfit=bs.mgfit, wor.c=wor.c, sp=sp) else{


    if((VC$fp==FALSE && (VC$l.sp1!=0 || VC$l.sp2!=0)) || PL!="P"){

       stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0  
       

	  while( stoprule.SP > pr.tolsp ){ 

             coefo <- fit$argument
             spEo <- spo <- sp 

             if(PL!="P" && (VC$l.sp1!=0 || VC$l.sp2!=0) ) spEo <- spo[-exclu]   
  
		 wor.c <- working.comp1(fit) 
  
                	bs.mgfit <- try(magic(y = wor.c$Z,  
                	                      X = wor.c$X,
                	                      sp= sp, S = qu.mag$Ss,
                        	              off = qu.mag$off, rank = qu.mag$rank,
                                	      gcv = FALSE,
                                	      gamma = VC$gamma))
                		if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                	spE <- sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1; names(sp) <- names(spo) 
                        if(PL!="P" && (VC$l.sp1!=0 || VC$l.sp2!=0) ) spE <- sp[-exclu]
                           
             o.ests <- c(fit$argument)

             fit <- trust(func.opt, o.ests, rinit=rinit, rmax=rmax,    
                              sp.xi1=sp["xi1"], sp.xi2=sp["xi2"], 
                              PL=PL, eqPL=eqPL, valPL=valPL, fitPL=fitPL, 
                              respvec=respvec, VC=VC, 
                              sp=spE, qu.mag=qu.mag, 
                              blather=TRUE, iterlim=iterlim)
                              
             iter.inner <- iter.inner + fit$iterations   	                                    
              

             n.ests <- c(fit$argument)

             if(iter.sp >= iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- max(abs(o.ests-n.ests))
                  
           
          }
          

       if(VC$gc.l == TRUE) gc()   
       
       magpp <- magic.post.proc(wor.c$X, bs.mgfit)
       

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


}





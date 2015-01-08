SemiParBIVProbit.fit  <- function(func.opt, start.v, 
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

  fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, parscale=VC$parscale,  
                sp.xi1=sp["xi1"], sp.xi2=sp["xi2"], PL=PL, eqPL=eqPL, 
                valPL=valPL, fitPL=fitPL, respvec=respvec, VC=VC,
                sp=spE, qu.mag=qu.mag, 
                blather=TRUE, iterlim = iterlim)  

  iter.if <- fit$iterations
  #n <- sum(as.numeric(fit$good)) 

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL

if( PL!="P" && fitPL!="pLiksp" ) list(fit=fit, iter.if=iter.if, conv.sp=conv.sp, iter.sp=iter.sp, iter.inner=iter.inner, bs.mgfit=bs.mgfit, wor.c=wor.c, sp=sp) else{



    if((VC$fp==FALSE && (VC$l.sp1!=0 || VC$l.sp2!=0)) || PL!="P"){

       stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0  
       
         if(PL == "P"){ 
       
         if(VC$awlm == FALSE) myf <- function(x) list( rbind( c(x[1],x[2]), 
                                                              c(x[2],x[3])  ) )
                                          
         if(VC$awlm == TRUE)  myf <- function(x) list( rbind( c(x[1],x[2],x[4]), 
                                                              c(x[2],x[3],x[5]),
                                                              c(x[4],x[5],x[6])  ) )                                                                
                      }
       
         if(PL != "P"){
       
              if(eqPL!="both")  myf <- function(x) list( rbind( c(x[1],x[2],x[4]),
                                                                c(x[2],x[3],x[5]),
                                                                c(x[4],x[5],x[6])  ) )
                    
              if(eqPL=="both")  myf <- function(x) list( rbind( c(x[1],x[2],x[4],x[5]),
                                                                c(x[2],x[3],x[6],x[7]),
                                                                c(x[4],x[6],x[8],x[9]),
                                                                c(x[5],x[7],x[9],x[10])  ) )
                      }
                      
              
              
if(PL == "P" && VC$awlm == TRUE){

  ll <- p.const <- NULL
  X  <- Matrix(0, 3*VC$n, (VC$X1.d2 + VC$X2.d2 + 1))
  D  <- Matrix(0, 3*VC$n, 1)
  
  wcf <- working.comp3

}else{              
              
  ll <- length(fit$argument) 
  p.const <- 2
  if(PL != "P"){   if(eqPL=="both"){ ll <- ll-2; p.const <- 4 } else {ll <- ll-1; p.const <- 3} } 
  
  if(PL == "P") X <- Matrix(0,p.const*VC$n,(VC$X1.d2+VC$X2.d2)) 
  if(PL != "P"){ if(eqPL == "both") X <- Matrix(0,p.const*VC$n,(VC$X1.d2+VC$X2.d2+2)) else X <- Matrix(0,p.const*VC$n,(VC$X1.d2+VC$X2.d2+1))   } 
  
  D <- Matrix(0, p.const*VC$n, 1) 
  
  wcf <- working.comp                      

}

  extra.l <- list(ll = ll, X = X, D = D, p.const = p.const)   
  rm(ll, X, D, p.const)
  if(VC$gc.l == TRUE) gc()

  

	  while( stoprule.SP > pr.tolsp ){ 

             coefo <- fit$argument
             spEo <- spo <- sp 
             if(PL!="P" && (VC$l.sp1!=0 || VC$l.sp2!=0) ) spEo <- spo[-exclu]   
  
		 wor.c <- try( wcf(fit, VC, myf, extra.l) ) 
		 

                 if(class(wor.c)=="try-error") break
             
                	bs.mgfit <- try(magic(y = wor.c$rW.Z, 
                	                      X = wor.c$rW.X,
                	                      sp= sp, S = qu.mag$Ss,
                        	              off = qu.mag$off, rank = qu.mag$rank,
                                	      gcv = FALSE,
                                	      gamma = VC$infl.fac))
                		if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                	spE <- sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1; names(sp) <- names(spo) 
                        if(PL!="P" && (VC$l.sp1!=0 || VC$l.sp2!=0) ) spE <- sp[-exclu]
                           
             o.ests <- c(fit$argument)

             fit <- try(trust(func.opt, o.ests, rinit=rinit, rmax=rmax, parscale=VC$parscale,
                              sp.xi1=sp["xi1"], sp.xi2=sp["xi2"], 
                              PL=PL, eqPL=eqPL, valPL=valPL, fitPL=fitPL, 
                              respvec=respvec, VC=VC, 
                              sp=spE, qu.mag=qu.mag, 
                              blather=TRUE, iterlim=iterlim),silent=TRUE) 

              			if(class(fit)=="try-error"){ 
              	 			fit  <- trust(func.opt, coefo, rinit=rinit, rmax=rmax, parscale=VC$parscale,
              	 			              sp.xi1=spo["xi1"], sp.xi2=spo["xi2"], 
              	 			              PL=PL, eqPL=eqPL, fitPL=fitPL,   
                                                      respvec=respvec, VC=VC,
                                                      sp=spEo, 
                		              	      qu.mag=qu.mag, blather=TRUE, iterlim=iterlim)
               		        conv.sp <- FALSE; break
                	                                    } 
             iter.inner <- iter.inner + fit$iterations   	                                    
              

             n.ests <- c(fit$argument)

             if(iter.sp>iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- max(abs(o.ests-n.ests))
                  
           
          }
          
       rm(extra.l)
       if(VC$gc.l == TRUE){ rm(wor.c); gc(); wor.c <- NULL}
       
       magpp <- magic.post.proc(wor.c$rW.X, bs.mgfit)
      
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





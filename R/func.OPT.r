func.OPT <- function(margins, M, type = "copR"){
  
  
if(type == "biv"){

  if(M$Model=="B" && margins[2] %in% M$bl  &&  is.null(M$theta.fx)) func.opt <- bprobgHs   
  if(M$Model=="B" && margins[2] %in% M$bl  && !is.null(M$theta.fx)) func.opt <- bprobgHsFixTheta

  if(M$Model=="B" && margins[2] %in% M$m1d  ) {func.opt <- bprobgHsDiscr1 ; func.optUniv <- bprobgHsContUniv}   
  if(M$Model=="B" && margins[2] %in% M$m2d  ) {func.opt <- bprobgHsDiscr2 ; func.optUniv <- bprobgHsContUniv} 
  
  if(M$Model=="BPO")                        func.opt <- bprobgHsPO 
  if(M$Model=="BPO0")                       func.opt <- bprobgHsPO0   
  if(M$Model=="BSS")                        func.opt <- bprobgHsSS   
  
  if(M$Model=="B" && margins[2] %in% M$m2  ) {func.opt <- bprobgHsCont  ; func.optUniv <- bprobgHsContUniv}   
  if(M$Model=="B" && margins[2] %in% M$m3  ) {func.opt <- bprobgHsCont3 ; func.optUniv <- bprobgHsContUniv3}  
  
}
  
  
  
if(type == "copR"){

  if(margins[1] %in% M$m1d && margins[2] %in% M$m2) func.opt  <- bdiscrcont12
  if(margins[1] %in% M$m1d && margins[2] %in% M$m3) func.opt  <- bdiscrcont13
  if(margins[1] %in% M$m2d && margins[2] %in% M$m2) func.opt  <- bdiscrcont
  if(margins[1] %in% M$m2d && margins[2] %in% M$m3) func.opt  <- bdiscrcont23
  
  if(margins[1] %in% M$m1d && margins[2] %in% M$m1d) func.opt  <- bdiscrdiscr11 
  if(margins[1] %in% M$m1d && margins[2] %in% M$m2d) func.opt  <- bdiscrdiscr12
  if(margins[1] %in% M$m2d && margins[2] %in% M$m2d) func.opt  <- bdiscrdiscr  
  
  
  if(margins[1] %in% M$m2 && margins[2] %in% M$m2 && M$BivD != "T") {if(M$robust == FALSE) func.opt <- bcont else func.opt <- bcontROB } 
  if(margins[1] %in% M$m2 && margins[2] %in% M$m2 && M$BivD == "T") func.opt <- bconttwoParC  

  if(margins[1] %in% M$m3 && margins[2] %in% M$m3 && M$BivD != "T") func.opt <- bcont3       
  if(margins[1] %in% M$m3 && margins[2] %in% M$m3 && M$BivD == "T") func.opt <- bcont3twoParC   
  
  if(margins[1] %in% M$m2 && margins[2] %in% M$m3 && M$BivD != "T") func.opt <- bcont23
  if(margins[1] %in% M$m2 && margins[2] %in% M$m3 && M$BivD == "T") func.opt <- bcont23twoParC
  
  if(margins[1] %in% M$m3 && margins[2] %in% M$m2 && M$BivD != "T") func.opt <- bcont32 
  if(margins[1] %in% M$m3 && margins[2] %in% M$m2 && M$BivD == "T") func.opt <- bcont32twoParC 


}


if(type == "copSS"){


  if(margins[2] %in% M$m2  ) func.opt <- bprobgHsContSS    
  if(margins[2] %in% M$m3  ) func.opt <- bprobgHsCont3SS  
  
  if(margins[2] %in% M$m1d  ) func.opt <- bprobgHsDiscr1SS
  if(margins[2] %in% M$m2d  ) func.opt <- bprobgHsDiscr2SS  


}
  

func.opt


}


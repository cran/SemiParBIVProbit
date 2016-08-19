eta.tr <- function(vrb.st, margin){

mupos <- c("LN","WEI","iG","GA","DAGUM","SM","FISK","NBI","NBII","PO","ZTP","PIG")
mub   <- c("BE")

if(margin %in% mupos){
   vrb.st <- ifelse( vrb.st > 20,   20, vrb.st )  # it was 28
   vrb.st <- ifelse( vrb.st < -9, -9, vrb.st ) # it was -17
}

if(margin %in% mub){
   vrb.st <- ifelse( vrb.st >  16,  16, vrb.st )  
   vrb.st <- ifelse( vrb.st < -16, -16, vrb.st ) 
}


vrb.st  
 
} 




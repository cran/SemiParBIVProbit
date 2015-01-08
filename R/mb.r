mb <- function(treat, outc, Model, B = 100, sig.lev = 0.05){


if(length(table(treat))!=2 ) stop("The treatment variable must be binary.")
if(length(table(outc))!=2  ) stop("The outcome variable must be binary.")
if(length(treat) != length(outc)  ) stop("The treatment and outcome variables have different lengths.")
if(!(Model %in% c("B", "BSS")) || missing(Model)) stop("Error in parameter Model value. It should be one of: B or BSS.")

IV <- NULL

datato <- data.frame(treat, outc)

if(!is.null(IV)){

	LBs <- UBs <- NULL
	t.s <- table(IV)
	if(length(t.s) == length(IV)) stop("Procedure not designed to do anything meaningful with this IV variable.")
	cat.names <- dimnames( t.s )[[1]]
	
}

n <- length(treat)
mind <- c(1:n)
UBb <- LBb <- NA 


################
# Main functions
################

if(Model == "B"){

mbf <- function(dat){

tab1 <- table(dat$treat)
pT1  <- prop.table(tab1)[2] 
pT0  <- 1 - pT1 
tab2 <- table(dat$treat, dat$outc)

pY1cT1 <- prop.table(tab2,1)[4] 
pY1cT0 <- prop.table(tab2,1)[3] 

stuff  <- pY1cT1*pT1 - pY1cT0*pT0

UB <- (stuff + pT0)
LB <- (stuff - pT1)

avte <- (pY1cT1 - pY1cT0)

list(UB = UB, LB = LB, avte = avte)

}


}


################

if(Model == "BSS"){

mbf <- function(dat){

tab1 <- table(dat$treat)
pSel1  <- prop.table(tab1)[2] 
pSel0  <- 1 - pSel1 
tab2 <- table(dat$treat, dat$outc)

pY1cS1 <- prop.table(tab2,1)[4] 

stuff  <- pY1cS1*pSel1

UB <- (stuff + pSel0)
LB <-  stuff

avte <- pY1cS1

list(UB = UB, LB = LB, avte = avte)

}


}

################

Cn.fun <- function(Cn, UB, LB, sigma.UB, sigma.LB, alpha, n){

pnorm(Cn + sqrt(n)*( (UB - LB)/max(sigma.LB, sigma.UB) ) ) - pnorm(-Cn) - alpha
    
}

################

if(!is.null(IV)){

mbfIV <- function(t.s, treat, outc, IV, cat.names){

for( i in 1:length(t.s) ){

mbfres <- mbf(data.frame(treat, outc)[IV==cat.names[i],])

LBs[i] <- mbfres$LB
UBs[i] <- mbfres$UB

}

LB <- max(LBs)
UB <- min(UBs)

list(LB = LB, UB = UB)

}

}
################
################






if(is.null(IV)){

########
# for sd
########

for(i in 1:B){

bootind  <- sample( mind, size = n, replace = TRUE)
mbfresB  <- mbf(datato[bootind,]) 

UBb[i] <- mbfresB$UB
LBb[i] <- mbfresB$LB

}

sUBb <- sd(UBb)
sLBb <- sd(LBb)

##############

mbfres <- mbf(datato)

avte <- as.numeric(format(mbfres$avte*100, digits = 3, trim=TRUE))
WCB  <- as.numeric(format(c(mbfres$LB*100, mbfres$UB*100), digits = 3, trim=TRUE))

Cn <- abs(uniroot( Cn.fun, interval=c(-5,5), UB = mbfres$UB, LB = mbfres$LB, sigma.UB = sUBb, sigma.LB = sLBb, alpha = sig.lev, n = n)$root) 
#Cn <- abs(nleqslv( 0.5, Cn.fun, UB = mbfres$UB, LB = mbfres$LB, sigma.UB = sUBb, sigma.LB = sLBb, alpha = sig.lev, n = n)$x) 


CIb1 <- as.numeric(round((mbfres$LB - Cn*sLBb/sqrt(n))*100, 2) )   
CIb2 <- as.numeric(round((mbfres$UB + Cn*sUBb/sqrt(n))*100, 2) )

out <- list(LB = WCB[1], UB = WCB[2], CI = c(CIb1,CIb2), av.p = avte, Model = Model) 

}













if(!is.null(IV)){

rei <- mbfIV(t.s, treat, outc, IV, cat.names)

LB <- as.numeric(format(rei$LB*100, digits = 3, trim=TRUE))
UB <- as.numeric(format(rei$UB*100, digits = 3, trim=TRUE))

avte <- as.numeric(format(mbf(data.frame(treat, outc))$avte*100, digits = 3, trim=TRUE))

########
# for sd
########

for(i in 1:B){

bootind  <- sample( mind, size = n, replace = TRUE)
bootsamp <- datato[bootind,]

IVbi <- IV[bootind]
t.s <- table(IVbi)
cat.names <- dimnames( t.s )[[1]]

mbfresB  <- mbfIV(t.s, bootsamp$treat, bootsamp$outc, IVbi, cat.names)

UBb[i] <- mbfresB$UB
LBb[i] <- mbfresB$LB

}

sUBb <- sd(UBb)
sLBb <- sd(LBb)

###############

Cn <- abs(uniroot( Cn.fun, interval=c(-5,5), UB = UB, LB = LB, sigma.UB = sUBb, sigma.LB = sLBb, alpha = sig.lev, n = n)$root) 


CIb1 <- as.numeric(round((rei$LB - Cn*sLBb/sqrt(n))*100, 2) )   
CIb2 <- as.numeric(round((rei$UB + Cn*sUBb/sqrt(n))*100, 2) )

out <- list(LB = LB, UB = UB, CI = c(CIb1,CIb2), av.p = avte, Model = Model) 



}


class(out) <- "mb"

out

}

